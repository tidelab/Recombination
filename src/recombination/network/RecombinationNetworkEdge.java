package recombination.network;

import java.util.*;

import beast.evolution.tree.Tree;

public class RecombinationNetworkEdge {

    public RecombinationNetworkNode parentNode, childNode;
    public BreakPoints breakPoints;
    public BreakPoints passingRange;
    public BreakPoints carryingRange;
    
    /**
     * status of this node after an operation is performed on the state *
     */
    int isDirty = Tree.IS_CLEAN;

    public boolean visited;
   
    // keeps track of the matrices for the likelihood calculations
    public double[] matrixList;

    public RecombinationNetworkEdge() { }

    public RecombinationNetworkEdge(RecombinationNetworkNode parentNode, RecombinationNetworkNode childNode,
    		BreakPoints breakPoints) {
        this.parentNode = parentNode;
        this.childNode = childNode;
        this.breakPoints = breakPoints;
    }

    public RecombinationNetworkEdge(RecombinationNetworkNode parentNode, RecombinationNetworkNode childNode,
    		BreakPoints breakPoints, double[] matrixList) {
        this.parentNode = parentNode;
        this.childNode = childNode;
        this.breakPoints = breakPoints;
        if (matrixList!=null) {
        	this.matrixList = new double[matrixList.length];
        	System.arraycopy(matrixList, 0, this.matrixList, 0, matrixList.length);
        }
        isDirty = Tree.IS_FILTHY;
    }
       
    public RecombinationNetworkEdge(RecombinationNetworkNode parentNode, RecombinationNetworkNode childNode,
    		int totalLength) {
        this.parentNode = parentNode;
        this.childNode = childNode;
        this.breakPoints = new BreakPoints(totalLength);
        isDirty = Tree.IS_FILTHY;
   }
    
    public RecombinationNetworkEdge(RecombinationNetworkNode parentNode, RecombinationNetworkNode childNode,
    		List<Integer> breakPointsList) {
        this.parentNode = parentNode;
        this.childNode = childNode;
        this.breakPoints = new BreakPoints();
        this.breakPoints.init(breakPointsList);
        isDirty = Tree.IS_FILTHY;
    }

    public double getRecombinationLength() {
        // There are always two reassortment configurations that
        // produce an unobserved reassortment: 1111 and 0000
        // (assuming 4 segs on lineage)
        return breakPoints.getLength();
    }

    public double getLength() {
        return parentNode.getHeight() - childNode.getHeight();
    }

    public boolean isRootEdge() {
        return parentNode == null;
    }
    
    public boolean isLeafEdge() {
    	return childNode.isLeaf();
    }

    public RecombinationNetworkEdge getCopy() {
        return getCopy(new HashMap<>());
    }

    public RecombinationNetworkEdge getCopy(Map<RecombinationNetworkNode,RecombinationNetworkNode> seenNodes) {
        RecombinationNetworkEdge edgeCopy;
       	edgeCopy = new RecombinationNetworkEdge(null, null, breakPoints.copy());
       	if (matrixList!=null) {
	       	edgeCopy.matrixList = new double[matrixList.length];
	    	System.arraycopy(matrixList, 0,  edgeCopy.matrixList, 0, matrixList.length);
       	}
       
        RecombinationNetworkNode childNodeCopy;
        boolean traverse = true;
        if (seenNodes.containsKey(childNode)) {
            childNodeCopy = seenNodes.get(childNode);
            traverse = false;
        } else {
            childNodeCopy = new RecombinationNetworkNode();
            childNodeCopy.setHeight(childNode.getHeight());
            childNodeCopy.setTaxonLabel(childNode.getTaxonLabel());
            childNodeCopy.setTaxonIndex(childNode.getTaxonIndex());
            childNodeCopy.setTypeIndex(childNode.typeIndex);
            childNodeCopy.setTypeLabel(childNode.typeLabel);
            childNodeCopy.setTypeLabel(childNode.typeLabel);
            if (childNode.states!=null) {
            	childNodeCopy.states = new int[childNode.states.length];
            	System.arraycopy(childNode.states, 0,  childNodeCopy.states, 0, childNode.states.length);
            }
            if (childNode.partials!=null) {
            	childNodeCopy.partials = new double[childNode.partials.length];
            	System.arraycopy(childNode.partials, 0,  childNodeCopy.partials, 0, childNode.partials.length);
            }
            seenNodes.put(childNode, childNodeCopy);
        }

        childNodeCopy.addParentEdge(edgeCopy);
       	edgeCopy.isDirty = isDirty;

        if (traverse) {
            for (RecombinationNetworkEdge childEdge : childNode.getChildEdges()) {
                RecombinationNetworkEdge childEdgeCopy = childEdge.getCopy(seenNodes);
                childNodeCopy.addChildEdge(childEdgeCopy);
            }
        }

        return edgeCopy;
    }
    
    /**
     * set the range of loci that goes left a this recombination node, only if 
     * @param from
     * @param to
     */
    public void setPassingRange(int from, int to) {
    	passingRange = new BreakPoints(from, to);
    }
    
    public BreakPoints getPassingRange() {
    	return passingRange;
    }

	public void setPassingRange(BreakPoints lociToDivert) {
		passingRange = lociToDivert;
	}

	
	public void setMatrix(double[] matrixList) {
		// TODO, allow for multiple
		System.arraycopy(this.matrixList, 0, matrixList,
                0, matrixList.length);
	}	
	
	public double[] getMatrix() {
		return matrixList;
	}
	
    public int isDirty() {
        return isDirty;
    }

    public void makeDirty(final int dirty) {
        isDirty |= dirty;
    }



	
}
