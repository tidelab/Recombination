package recombination.operators;

import beast.core.Input;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import recombination.network.BreakPoints;
import recombination.network.RecombinationNetworkEdge;
import recombination.network.RecombinationNetworkNode;

import java.util.BitSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class DivertLociOperator extends EmptyEdgesRecombinationNetworkOperator {

    public Input<Double> scaleFactorInput = new Input<>(
            "scaleFactor",
            "Scale factor tuning parameter.",
            1.0);

    double lambdaDiversion;
    public int totalLength;
    
    boolean stop = false;
    
    public void initAndValidate() {
    	lambdaDiversion = scaleFactorInput.get();    	
    	super.initAndValidate();
    	totalLength = network.totalLength;
    }
	
    @Override
    public double networkProposal() {    	    	
        double logHR = 0.0;
        
        List<RecombinationNetworkNode> sourceNodes = network.getNodes().stream()
                .filter(e -> e.isRecombination())
                .filter(e -> !e.getParentEdges().get(0).breakPoints.isEmpty())
                .filter(e -> !e.getParentEdges().get(1).breakPoints.isEmpty())
                .collect(Collectors.toList());

        if (sourceNodes.isEmpty())
            return Double.NEGATIVE_INFINITY;

        logHR -= Math.log(1.0/sourceNodes.size());

        network.startEditing(this);
        
    	int newBreakPoint = Randomizer.nextInt(totalLength-1)+1;
       
        // get both edges
        RecombinationNetworkNode node = sourceNodes.get(Randomizer.nextInt(sourceNodes.size()));
        RecombinationNetworkEdge edge1,edge2;
        if (node.getParentEdges().get(0).breakPoints.getMax()>node.getParentEdges().get(1).breakPoints.getMax()) {
        	edge1 = node.getParentEdges().get(1);
        	edge2 = node.getParentEdges().get(0);
        }else {
        	edge1 = node.getParentEdges().get(0);
        	edge2 = node.getParentEdges().get(1);       	
        }
              

        edge1.makeDirty(Tree.IS_DIRTY);
        edge2.makeDirty(Tree.IS_DIRTY);

    	edge2.passingRange = new BreakPoints(newBreakPoint, totalLength-1);
    	edge1.passingRange = new BreakPoints(0,newBreakPoint-1);

    	if (edge2.breakPoints.getMin()>newBreakPoint) {
        	BreakPoints rangeToDivert = new BreakPoints(newBreakPoint, edge2.breakPoints.getMin()-1);
        	rangeToDivert.and(edge1.breakPoints);
	        logHR -= addLociToAncestors(edge2, rangeToDivert);
	        logHR += removeLociFromAncestors(edge1, rangeToDivert);
        }else {
        	BreakPoints rangeToDivert = new BreakPoints(edge2.breakPoints.getMin(), newBreakPoint-1);
        	rangeToDivert.and(edge2.breakPoints);
	        logHR -= addLociToAncestors(edge1, rangeToDivert);
	        logHR += removeLociFromAncestors(edge2, rangeToDivert);
        }
        
        int reverseSourceEdgeCount = (int)(network.getNodes().stream()
                .filter(e -> e.isRecombination())
                .filter(e -> !e.getParentEdges().get(0).breakPoints.isEmpty())
                .filter(e -> !e.getParentEdges().get(1).breakPoints.isEmpty())
                .count());

        logHR += Math.log(1.0/reverseSourceEdgeCount);
        return logHR;
    }


    /**
     * Remove segments from this edge and ancestors.
     *
     * @param edge edge at which to start removal
     * @param segsToRemove segments to remove from edge and ancestors
     * @return log probability of reverse operation
     */
    public double removeLociFromAncestors(RecombinationNetworkEdge edge, BreakPoints rangeToRemove) {
        double logP = 0.0;
        

        rangeToRemove = rangeToRemove.copy();
        
        rangeToRemove.and(edge.breakPoints);
        
        if (rangeToRemove.isEmpty())
            return logP;
        
        if (edge.isRootEdge())
            return logP;

        edge.childNode.dirtyBreakPoints = new BreakPoints(0,totalLength-1);

        edge.breakPoints.andNot(rangeToRemove);                      
        
        edge.makeDirty(Tree.IS_FILTHY); 

        if (edge.parentNode.isRecombination()) {
        	
//        	// get the breakpoints before the operation
//        	int old_max1 = edge.parentNode.getParentEdges().get(0).breakPoints.getMax();
//        	int old_max2 = edge.parentNode.getParentEdges().get(1).breakPoints.getMax();
        	
            logP += removeLociFromAncestors(edge.parentNode.getParentEdges().get(0), rangeToRemove);
            logP += removeLociFromAncestors(edge.parentNode.getParentEdges().get(1), rangeToRemove);
            
//    		int min1 = edge.parentNode.getParentEdges().get(0).breakPoints.getMin();
//    		int max1 = edge.parentNode.getParentEdges().get(0).breakPoints.getMax();
//    		
//    		int min2 = edge.parentNode.getParentEdges().get(1).breakPoints.getMin();
//    		int max2 = edge.parentNode.getParentEdges().get(1).breakPoints.getMax();
//    		
//            if (min1==-1 && min2==-1) {
//            	// both passing ranges are null
//            	logP += Math.log(0.5) + Math.log(1.0/(totalLength));
////            	return Double.NEGATIVE_INFINITY;
//            }else if (min1==-1){ //corr
//            	logP += getProbNewPassingRange(min2, max2, old_max1, rangeToRemove);
////            	return Double.NEGATIVE_INFINITY;
//            }else if(min2==-1) { //corr
//            	logP += getProbNewPassingRange(min1, max1, old_max2, rangeToRemove);
////            	return Double.NEGATIVE_INFINITY;
//            }else {
//            	// check if the parts of the removed breakpoints are between the breakpoints of the two edges
////            	BreakPoints bp = new BreakPoints(Math.min(max1, max2), Math.max(min1, min2));
////            	int diff = Math.max(min1, min2) - Math.min(max1, max2);
////            	logP += Math.log(1.0/(diff));            	     	
//            }  
        } else {
        	rangeToRemove.andNot(getSisterEdge(edge).breakPoints);
            logP += removeLociFromAncestors(edge.parentNode.getParentEdges().get(0), rangeToRemove);
        }

        return logP;
    }

    /**
     * Add segments to this edge and ancestors.
     *
     * @param edge edge at which to start addition
     * @param segsToAdd segments to add to the edge and ancestors
     * @return log probability of operation
     */
    public double addLociToAncestors(RecombinationNetworkEdge edge, BreakPoints rangeToAdd) {
        double logP = 0.0;

        rangeToAdd = rangeToAdd.copy();
               
        if (rangeToAdd.isEmpty())
            return logP;

        rangeToAdd.andNot(edge.breakPoints);

        if (rangeToAdd.isEmpty())
            return logP;        


        edge.breakPoints.or(rangeToAdd);
        
        edge.makeDirty(Tree.IS_FILTHY); 


        if (edge.isRootEdge())
            return logP;        
        
        if (edge.parentNode.isRecombination()) {        	
        	// resample the passing Range between the boundries given by the left and right breakpoints
//        	logP += resamplePassingRange(edge.parentNode.getParentEdges().get(0),
//        			edge.parentNode.getParentEdges().get(1), rangeToAdd);  

            BreakPoints rangeToAddLeft = rangeToAdd.copy();
            BreakPoints rangeToAddRight = rangeToAdd.copy();
                        
            rangeToAddLeft.and(edge.parentNode.getParentEdges().get(0).passingRange);
            rangeToAddRight.and(edge.parentNode.getParentEdges().get(1).passingRange);
                        
            logP += addLociToAncestors(edge.parentNode.getParentEdges().get(0), rangeToAddLeft);
            logP += addLociToAncestors(edge.parentNode.getParentEdges().get(1), rangeToAddRight);
        } else {
            logP += addLociToAncestors(edge.parentNode.getParentEdges().get(0), rangeToAdd);
        }
        return logP;
    }
       
       
    private double resamplePassingRange(RecombinationNetworkEdge edge1,
			RecombinationNetworkEdge edge2, BreakPoints rangeToAdd) {
    	
    	double logHR = 0.0;

		int min1 = edge1.breakPoints.getMin();
		int max1 = edge1.breakPoints.getMax();
		
		int min2 = edge2.breakPoints.getMin();
		int max2 = edge2.breakPoints.getMax();	

		if (min1==-1 && min2==-1) { 
			// both passing ranges are null
			if (Randomizer.nextBoolean()) {
				sampleNullRange(edge1, edge2);
			}else {
				sampleNullRange(edge2, edge1);
			}	
			logHR += Math.log(0.5) + Math.log(1.0/(totalLength));
		}else if(min1==-1) {
			// 1 is null
			logHR += sampleNewPassingRange(edge2,edge1,min2,max2);
		}else if(min2==-1) {
			// 2 is null
			logHR += sampleNewPassingRange(edge1,edge2,min1,max1);
		}else {		
//			// resample between the two edges
//			if (max1 > max2) {
//				logHR += sampleBetweenRange(edge1, edge2, min1, max2, rangeToAdd);
//			}else {
//				logHR += sampleBetweenRange(edge2, edge1, min2, max1, rangeToAdd);
//			}
		}   
		return logHR;
	}
    
    private void sampleNullRange(RecombinationNetworkEdge edge1, RecombinationNetworkEdge edge2) {
		int start = Randomizer.nextInt(totalLength);
		if (start==0) {
			edge1.setPassingRange(start, totalLength-1);
			edge2.setPassingRange(null);
		}else {
			edge1.setPassingRange(start, totalLength-1);
			edge2.setPassingRange(0,start-1);
		}    	
    }
    
    private double sampleBetweenRange(RecombinationNetworkEdge edge1, RecombinationNetworkEdge edge2, 
    		int min1, int max2, BreakPoints rangeToAdd) {
    	int diff = min1 - max2;
		int newBreakPoint = Randomizer.nextInt(diff)+max2;
	
		edge1.setPassingRange(newBreakPoint+1, totalLength-1);
		edge2.setPassingRange(0, newBreakPoint);
    	return Math.log(1.0/(diff));
    }
        
    private double sampleNewPassingRange(RecombinationNetworkEdge edge1,
			RecombinationNetworkEdge edge2, int min, int max) {
    	double logHR = 0.0;

    	
		if (min==0 && max==totalLength-1) {
			// passing range 2 is null
			edge1.setPassingRange(0, totalLength-1);
			edge2.passingRange=null;
		}else if(min==0) {
			
			int diff = totalLength-max;
			logHR += Math.log(1.0/diff);				
			int start = Randomizer.nextInt(diff)+max;
				
			if (start==(totalLength-1)) 
				edge2.passingRange = null;
			else
				edge2.setPassingRange(start+1, totalLength-1);
			
			edge1.setPassingRange(0,start);
		}else if(max==totalLength-1) {
			logHR += Math.log(1.0/(min+1));

			int end = Randomizer.nextInt(min+1)-1;
			if (end==-1) 
				edge2.setPassingRange(null);
			else 
				edge2.setPassingRange(0, end);			

			edge1.setPassingRange(end+1,totalLength-1);						
    	}else {
    		logHR += Math.log(0.5);
			if (Randomizer.nextBoolean()) {		
				int diff = totalLength-max;
				logHR += Math.log(1.0/diff);				
				
				int start = Randomizer.nextInt(diff)+max;
				
				if (start==(totalLength-1)) 
					edge2.passingRange = null;
				else
					edge2.setPassingRange(start+1, totalLength-1);
				
				edge1.setPassingRange(0,start);
			}else {
				logHR += Math.log(1.0/(min+1));
				int end = Randomizer.nextInt(min+1)-1;			
				
				if (end==-1)
					edge2.passingRange = null;
				else
					edge2.setPassingRange(0, end);
				
				edge1.setPassingRange(end+1,totalLength-1);
			}	
		
		}   

    	return logHR;
    }

    private double getProbNewPassingRange(double min, double max, double old_max , BreakPoints bp) {

		if (min==0 && max==totalLength-1) {
			return 0.0;
    	}else if (min==0) {
    		return Math.log(1.0/(totalLength-max));
		}else if (max==totalLength-1) {
			return Math.log(1.0/(min+1));
    	}else {
        	if (old_max>=max) {
        		return Math.log(0.5) + Math.log(1.0/(totalLength-max));
        	}else if (old_max<min) {
        		return Math.log(0.5) + Math.log(1.0/(min+1));
        	}else {
        		throw new IllegalArgumentException("should not happen, or at least is not accounted for");
        	}
    	}
    }
}
