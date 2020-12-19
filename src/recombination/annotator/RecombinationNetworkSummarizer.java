/*
 * Copyright (C) 2015 Tim Vaughan <tgvaughan@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package recombination.annotator;

import beast.core.util.Log;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetwork;
import recombination.network.RecombinationNetworkNode;
import recombination.statistics.DotConverter;

import javax.swing.*;
import javax.swing.border.EtchedBorder;
import java.awt.*;
import java.io.*;
import java.lang.reflect.InvocationTargetException;
import java.util.*;
import java.util.List;

/**
 * A rewrite of TreeAnnotator targeted at summarizing ACG logs
 * generated by bacter.
 * 
 * @author Tim Vaughan <tgvaughan@gmail.com>
 * @author Nicola Felix Müller <nicola.felix.mueller@gmail.com>
 */
public class RecombinationNetworkSummarizer extends RecombinationAnnotator {

    private enum SummaryStrategy { MEAN, MEDIAN }

    private static class NetworkAnnotatorOptions {
        File inFile;
        File outFile = new File("summary.tree");
        File targetFile;
        double burninPercentage = 10.0;
        SummaryStrategy summaryStrategy = SummaryStrategy.MEAN;
        BreakPoints breakPoints = new BreakPoints();
        boolean useDotFormat = false;
        int useEveryTree = 0;

        @Override
        public String toString() {
            return "Active options:\n" +
                    "Input file: " + inFile + "\n" +
                    "Output file: " + outFile + "\n" +
                    "Target Network file: " + targetFile + "\n" +
                    "Burn-in percentage: " + burninPercentage + "\n" +
                    "Node height and conv. site summary: " + summaryStrategy + "\n" +
            		"Remove Loci for summary: " + breakPoints + "\n" +
            		"dot format output: " + useDotFormat + "\n" + 
    				"treesFormat output: " + useEveryTree + "\n";
       }
    }

    public RecombinationNetworkSummarizer(NetworkAnnotatorOptions options) throws IOException {

        // Display options:
        System.out.println(options + "\n");

        // Initialise reader

        RecombinationLogReader logReader = new RecombinationLogReader(options.inFile,
                options.burninPercentage);

        System.out.println(logReader.getNetworkCount() + " Networks in file.");

        System.out.println("The first " + logReader.getBurnin() +
                 " (" + options.burninPercentage + "%) ACGs will be discarded " +
                "to account for burnin.");
        
        // keeps track of the mcc network
        RecombinationNetwork bestNetwork = null;
        // keeps track of the leave nodes
        List<RecombinationNetworkNode> leafNodes = new ArrayList<>();
        
        boolean onTarget = options.targetFile == null ? false : true;
        
        if (!onTarget){
	        // get the clades for each reassortment event in every network
	       RecombinationNetworkCladeSystem cladeSystem = new RecombinationNetworkCladeSystem();
	        
	        // build the clades
	        boolean first = true;
	        for (RecombinationNetwork network : logReader){
	        	if (first){        		
	            	for (RecombinationNetworkNode networkNode : network.getNodes()){
	            		if (networkNode.isLeaf()){
	            			leafNodes.add(networkNode);
	            		}
	        		}
	            	cladeSystem.setLeafLabels(leafNodes, network.totalLength);
	        		first = false;
	        	}
	        	pruneNetwork(network, options.breakPoints);
	        	cladeSystem.add(network, true); 
	        }
	        
	        System.out.println("\nComputing CF clade credibilities...");
	        // calculate the network clade credibilities      
	        
	        // get the network with the highest count
	        double bestScore = Double.NEGATIVE_INFINITY;
	
	        for (RecombinationNetwork network : logReader ) {
	        	pruneNetwork(network, options.breakPoints);
	        	double score = cladeSystem.getLogCladeCredibility(network, logReader.getCorrectedNetworkCount());
	        	if (score>bestScore) {
	        		bestNetwork = network;
	        		bestScore = score;
	        	}
	        }
        }else{
	        System.out.println("\nRead in target network...");

            RecombinationLogReader targetLogReader = new RecombinationLogReader(options.targetFile, 0.0);
            int c= 0;
	        for (RecombinationNetwork network : targetLogReader ) {
	        	pruneNetwork(network, options.breakPoints);
	        	bestNetwork = network;
	        	c++;
	        }
	        if (c!=1)
	        	throw new IllegalArgumentException("more than one network found as target network");
	        
        	for (RecombinationNetworkNode networkNode : bestNetwork.getNodes()){
        		if (networkNode.isLeaf()){
        			leafNodes.add(networkNode);
        		}
    		}
        }
        

        // get the posterior probabilities of each coalescent network node
        RecombinationNetworkCladeSystem bestCladeSystem = new RecombinationNetworkCladeSystem();
//                
        // add leafnodes
    	bestCladeSystem.setLeafLabels(leafNodes, bestNetwork.totalLength);
    	    	
    	// build clade system
        bestCladeSystem.add(bestNetwork, true);
//        
        Set<String> attributeNames = new HashSet<>();
		attributeNames.add("height");
//		
        // print the network to file
        System.out.println("\nCollect Atributes...");
        for (RecombinationNetwork network : logReader ) {
        	pruneNetwork(network, options.breakPoints);
    		bestCladeSystem.collectAttributes(network, attributeNames, true);
    	}
        
        // print the network to file
        System.out.println("\nSummarize Atributes...");        
    	if (options.summaryStrategy == SummaryStrategy.MEAN)
    		bestCladeSystem.summarizeAttributes(bestNetwork, attributeNames, true, logReader.getCorrectedNetworkCount(), onTarget);
    	else
    		bestCladeSystem.summarizeAttributes(bestNetwork, attributeNames, false, logReader.getCorrectedNetworkCount(), onTarget);
//
    	// print the network to file
        System.out.println("\nWriting output to " + options.outFile.getName()
        	+ "...");
        try (PrintStream ps = new PrintStream(options.outFile)) {
        	if (options.useDotFormat) {
        		 List<String> dotString = DotConverter.getDotFormat(bestNetwork);
        		 for (String s : dotString) {
     	        	ps.print(s);
        		 }
        	}else if (options.useEveryTree>0) {
	        	ps.print(logReader.getPreamble());
	        	for (int i = 0; i < bestNetwork.totalLength; i=i+options.useEveryTree) {
	        		RecombinationNetwork n = new RecombinationNetwork(bestNetwork.getRootEdge().getCopy());
	        		n.totalLength = bestNetwork.totalLength;
	        		BreakPoints bp2 = new BreakPoints(i,i);
		        	pruneNetwork(n, bp2);
		        	
		        	
	        		ps.println("tree Loci_" + i + " = " + n.getExtendedNewickVerbose());
	        	}
	
	        	String postamble = logReader.getPostamble();
	        	if (postamble.length() > 0)
	        		ps.println(postamble);
	        	else
	        		ps.println("End;");

        	} else {
	        	ps.print(logReader.getPreamble());
	        	ps.println("tree STATE_0 = " + bestNetwork.getExtendedNewickVerbose());
	
	        	String postamble = logReader.getPostamble();
	        	if (postamble.length() > 0)
	        		ps.println(postamble);
	        	else
	        		ps.println("End;");
        	}
        }        

        System.out.println("\nDone!");
    }      
 
    /**
     * Use a GUI to retrieve ACGAnnotator options.
     *
     * @param options options object to populate using GUI
     * @return true if options successfully collected, false otherwise
     */
    private static boolean getOptionsGUI(NetworkAnnotatorOptions options) {

        boolean[] canceled = {false};

        JDialog dialog = new JDialog((JDialog)null, true);
        dialog.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        dialog.setLocationRelativeTo(null);
        dialog.setTitle("Reassortment Network Summarizer");

        JLabel logFileLabel = new JLabel("Reassortment Network log file:");
        JLabel outFileLabel = new JLabel("Output file:");
        JLabel targetFileLabel = new JLabel("Target file:");
        JLabel burninLabel = new JLabel("Burn-in percentage:");
        JLabel summaryMethodLabel = new JLabel("Position summary method:");
        JLabel removeSegmentLabel = new JLabel("Removes segments from the summary:");

        JTextField inFilename = new JTextField(20);
        inFilename.setEditable(false);
        JButton inFileButton = new JButton("Choose File");

        JTextField outFilename = new JTextField(20);
        outFilename.setText(options.outFile.getName());
        outFilename.setEditable(false);
        JButton outFileButton = new JButton("Choose File");
        
        JTextField targetFilename = new JTextField(20);
        targetFilename.setEditable(false);
        JButton targetFileButton = new JButton("Choose File");

        JSlider burninSlider = new JSlider(JSlider.HORIZONTAL,
                0, 100, (int)(options.burninPercentage));
        burninSlider.setMajorTickSpacing(50);
        burninSlider.setMinorTickSpacing(10);
        burninSlider.setPaintTicks(true);
        burninSlider.setPaintLabels(true);
        burninSlider.setSnapToTicks(true);
        
        JTextField removeSegments = new JTextField(20);
        removeSegments.setText("none");
        removeSegments.setEditable(true);


        JComboBox<SummaryStrategy> heightMethodCombo = new JComboBox<>(SummaryStrategy.values());

        Container cp = dialog.getContentPane();
        BoxLayout boxLayout = new BoxLayout(cp, BoxLayout.PAGE_AXIS);
        cp.setLayout(boxLayout);

        JPanel mainPanel = new JPanel();

        GroupLayout layout = new GroupLayout(mainPanel);
        mainPanel.setLayout(layout);
        layout.setAutoCreateGaps(true);
        layout.setAutoCreateContainerGaps(true);

        layout.setHorizontalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(logFileLabel)
                        .addComponent(outFileLabel)
                        .addComponent(targetFileLabel)
                        .addComponent(burninLabel)
                        .addComponent(summaryMethodLabel)
                        .addComponent(removeSegmentLabel))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFilename)
                        .addComponent(outFilename)
                        .addComponent(targetFilename)
                        .addComponent(burninSlider)
                        .addComponent(heightMethodCombo)
                        .addComponent(removeSegments))
                .addGroup(layout.createParallelGroup(GroupLayout.Alignment.LEADING, false)
                        .addComponent(inFileButton)
                        .addComponent(outFileButton)
                        .addComponent(targetFileButton)));

        layout.setVerticalGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup()
                        .addComponent(logFileLabel)
                        .addComponent(inFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(inFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(outFileLabel)
                        .addComponent(outFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(outFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(targetFileLabel)
                        .addComponent(targetFilename,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)
                        .addComponent(targetFileButton))
                .addGroup(layout.createParallelGroup()
                        .addComponent(burninLabel)
                        .addComponent(burninSlider,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE))
                .addGroup(layout.createParallelGroup()
                        .addComponent(summaryMethodLabel)
                        .addComponent(heightMethodCombo,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE))
                .addGroup(layout.createParallelGroup()
                        .addComponent(removeSegmentLabel)
                        .addComponent(removeSegments,
                                GroupLayout.PREFERRED_SIZE,
                                GroupLayout.DEFAULT_SIZE,
                                GroupLayout.PREFERRED_SIZE)));

        mainPanel.setBorder(new EtchedBorder());
        cp.add(mainPanel);

        JPanel buttonPanel = new JPanel();

        JButton runButton = new JButton("Analyze");
        runButton.addActionListener((e) -> {
            options.burninPercentage = burninSlider.getValue();
            options.summaryStrategy = (SummaryStrategy)heightMethodCombo.getSelectedItem();
//            if (!removeSegments.getText().contains("none")){
//            	String[] splitstr = removeSegments.getText().split("");
//            	options.removeSegments = new int[splitstr.length];
//            	for (int i = 0; i < splitstr.length; i++)
//            		options.removeSegments[i] = Integer.parseInt(splitstr[i]);
//            }
            dialog.setVisible(false);
        });
        runButton.setEnabled(false);
        buttonPanel.add(runButton);

        JButton cancelButton = new JButton("Quit");
        cancelButton.addActionListener((e) -> {
            dialog.setVisible(false);
            canceled[0] = true;
        });
        buttonPanel.add(cancelButton);

        JFileChooser inFileChooser = new JFileChooser();
        inFileButton.addActionListener(e -> {
            inFileChooser.setDialogTitle("Select ACG log file to summarize");
            if (options.inFile == null)
                inFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));
            int returnVal = inFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.inFile = inFileChooser.getSelectedFile();
                inFilename.setText(inFileChooser.getSelectedFile().getName());
                runButton.setEnabled(true);
            }
        });

        JFileChooser outFileChooser = new JFileChooser();
        outFileButton.addActionListener(e -> {
            outFileChooser.setDialogTitle("Select output file name.");
            if (options.inFile != null)
                outFileChooser.setCurrentDirectory(options.inFile);
            else
                outFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));

            outFileChooser.setSelectedFile(options.outFile);
            int returnVal = outFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.outFile = outFileChooser.getSelectedFile();
                outFilename.setText(outFileChooser.getSelectedFile().getName());
            }
        });
        
        JFileChooser targetFileChooser = new JFileChooser();
        targetFileButton.addActionListener(e -> {
            targetFileChooser.setDialogTitle("Select output file name.");
            if (options.inFile != null)
            	targetFileChooser.setCurrentDirectory(options.inFile);
            else
            	targetFileChooser.setCurrentDirectory(new File(System.getProperty("user.dir")));

            targetFileChooser.setSelectedFile(options.outFile);
            int returnVal = targetFileChooser.showOpenDialog(dialog);

            if (returnVal == JFileChooser.APPROVE_OPTION) {
                options.targetFile = targetFileChooser.getSelectedFile();
                outFilename.setText(targetFileChooser.getSelectedFile().getName());
            }
        });
             			
        cp.add(buttonPanel);

        dialog.pack();
        dialog.setResizable(false);
        dialog.setVisible(true);

        return !canceled[0];
    }

    /**
     * Prepare JFrame to which ACGAnnotator output streams will be
     * directed.
     */
    private static void setupGUIOutput() {

        JFrame frame = new JFrame();
        frame.setTitle("Reassortment Network Annotator");
        frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

        JTextArea textArea = new JTextArea(25, 80);
        textArea.setFont(new Font("monospaced", Font.PLAIN, 12));
        textArea.setEditable(false);
        frame.getContentPane().add(new JScrollPane(textArea), BorderLayout.CENTER);

        JButton closeButton = new JButton("Close");
        closeButton.addActionListener(e -> System.exit(0));
        JPanel buttonPanel = new JPanel();
        buttonPanel.add(closeButton);
        frame.getContentPane().add(buttonPanel, BorderLayout.PAGE_END);

        // Redirect streams to output window:
        OutputStream out = new OutputStream() {
            @Override
            public void write(int b) throws IOException {
                SwingUtilities.invokeLater(() -> {
                    if ((char)b == '\r') {
                        int from = textArea.getText().lastIndexOf("\n") + 1;
                        int to = textArea.getText().length();
                        textArea.replaceRange(null, from, to);
                    } else
                        textArea.append(String.valueOf((char) b));
                });
            }
        };

        System.setOut(new PrintStream(out, true));
        System.setErr(new PrintStream(out, true));

        frame.pack();
        frame.setVisible(true);
    }

    public static String helpMessage =
            "ACGAnnotator - produces summaries of Bacter ACG log files.\n"
                    + "\n"
                    + "Usage: appstore ACGAnnotator [-help | [options] logFile [outputFile]\n"
                    + "\n"
                    + "Option                   Description\n"
                    + "--------------------------------------------------------------\n"
                    + "-help                    Display usage info.\n"
                    + "-positions {mean,median} Choose position summary method.\n"
                    + "                         (default mean)\n"
                    + "-burnin percentage       Choose _percentage_ of log to discard\n"
                    + "                         in order to remove burn-in period.\n"
                    + "                         (Default 10%)\n"
                    + "-threshold percentage    Choose minimum posterior probability\n"
                    + "                         for including conversion in summary.\n"
                    + "                         (Default 50%)\n"
                    + "-recordGeneFlow gfFile   Record posterior distribution of gene\n"
                    + "                         flow in given file.\n"
                    + "\n"
                    + "If no output file is specified, output is written to a file\n"
                    + "named 'summary.tree'.";

    /**
     * Print usage info and exit.
     */
    public static void printUsageAndExit() {
        System.out.println(helpMessage);
        System.exit(0);
    }

    /**
     * Display error, print usage and exit with error.
     */
    public static void printUsageAndError(String errMsg) {
        System.err.println(errMsg);
        System.err.println(helpMessage);
        System.exit(1);
    }

    /**
     * Retrieve ACGAnnotator options from command line.
     *
     * @param args command line arguments
     * @param options object to populate with options
     */
    public static void getCLIOptions(String[] args, NetworkAnnotatorOptions options) {
        int i=0;
        while (args[i].startsWith("-")) {
            switch(args[i]) {
                case "-help":
                    printUsageAndExit();
                    break;

                case "-burnin":
                    if (args.length<=i+1)
                        printUsageAndError("-burnin must be followed by a number (percent)");

                    try {
                        options.burninPercentage = Double.parseDouble(args[i+1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("Error parsing burnin percentage.");
                    }

                    if (options.burninPercentage<0 || options.burninPercentage>100) {
                        printUsageAndError("Burnin percentage must be >= 0 and < 100.");
                    }

                    i += 1;
                    break;

                case "-positions":
                    if (args.length<=i+1) {
                        printUsageAndError("-positions must be followed by either 'MEAN' or 'MEDIAN'.");
                    }

                    if (args[i+1].toLowerCase().equals("mean")) {
                        options.summaryStrategy = SummaryStrategy.MEAN;

                        i += 1;
                        break;
                    }

                    if (args[i+1].toLowerCase().equals("median")) {
                        options.summaryStrategy = SummaryStrategy.MEDIAN;

                        i += 1;
                        break;
                    }

                    printUsageAndError("-positions must be followed by either 'MEAN' or 'MEDIAN'.");

                case "-subsetRange":
                    if (args.length<=i+1) {
                        printUsageAndError("-subsetRange must be a range in the format of 0-100.");
                    }

                    try {
                    	String[] argarray = args[i + 1].split(",");
                    	List<Integer> bp_list = new ArrayList<>();
                    	for (int j = 0; j < argarray.length; j++) {
                    		String[] tmp = argarray[j].split("-");
                    		bp_list.add(Integer.parseInt(tmp[0]));
                    		bp_list.add(Integer.parseInt(tmp[1]));
                    	}
                		options.breakPoints.init(bp_list);
                    } catch (NumberFormatException e) {
                        printUsageAndError("removeSegments must be an array of integers separated by commas if more than one");
                     }

                    i += 1;
                    break;
                    
                case "-target":
                    if (args.length<=i+1) {
                        printUsageAndError("-target must be followed by a network file.");
                    }

                    try {
                		options.targetFile = new File(args[i + 1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("removeSegments must be an array of integers separated by commas if more than one");
                     }

                    i += 1;
                    break;
                    
                case "-dotFormat":
                    if (args.length<=i+1) {
                        printUsageAndError("-dotFormat must be followed by true or false.");
                    }

                    try {
                		options.useDotFormat = Boolean.parseBoolean(args[i + 1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("dotFormat must be followed by true or false");
                     }

                    i += 1;
                    break;
                    
                case "-useEveryTree":
                    if (args.length<=i+1) {
                        printUsageAndError("-useEveryTree must be followed by true or false.");
                    }

                    try {
                		options.useEveryTree = Integer.parseInt(args[i + 1]);
                    } catch (NumberFormatException e) {
                        printUsageAndError("treesFormat must be followed by true or false");
                     }

                    i += 1;
                    break;


                default:
                    printUsageAndError("Unrecognised command line option '" + args[i] + "'.");
            }

            i += 1;
        }

        if (i >= args.length)
            printUsageAndError("No input file specified.");
        else
            options.inFile = new File(args[i]);

        if (i+1<args.length)
            options.outFile = new File(args[i+1]);
    }

    /**
     * Main method for ACGAnnotator.  Sets up GUI if needed then
     * uses the ACGAnnotator constructor to actually perform the analysis.
     *
     * @param args command line arguments
     */
    public static void main(String[] args) {
    	NetworkAnnotatorOptions options = new NetworkAnnotatorOptions();

        if (args.length == 0) {
            // Retrieve options from GUI:

            try {
                UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
            } catch (ClassNotFoundException | InstantiationException | UnsupportedLookAndFeelException | IllegalAccessException e) {
                Log.warning.println("Error setting cross-platform look and feel.");
            }

            try {
                SwingUtilities.invokeAndWait(() -> {
                    if (!getOptionsGUI(options))
                        System.exit(0);

                    setupGUIOutput();
                });
            } catch (InterruptedException | InvocationTargetException e) {
                e.printStackTrace();
            }
        } else {
            getCLIOptions(args, options);
        }

        // Run ACGAnnotator
        try {
            new RecombinationNetworkSummarizer(options);

        } catch (Exception e) {
            if (args.length == 0) {
                JOptionPane.showMessageDialog(null, e.getMessage(),
                        "Error", JOptionPane.ERROR_MESSAGE);
            } else {
                System.err.println("Error: " + e.getMessage());
                e.printStackTrace();
                System.err.println();
                System.err.println(helpMessage);
            }

            System.exit(1);
        }
    }
}