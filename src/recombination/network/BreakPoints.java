package recombination.network;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import beast.evolution.tree.Node;

/**
 * Keep track of recombination break points in the network
 * @author Nicola Müller
 */
public class BreakPoints {
	
	public List<Range> breakPoints;	
	private BreakPoints leftBreakPoints;
	private BreakPoints rightBreakPoints;    
	
	final static RangeComparator rc = new RangeComparator();

		
	public BreakPoints(int totalLength) { 
		breakPoints = new ArrayList<>();
		breakPoints.add(new Range(0, totalLength-1));
	}

	
	public BreakPoints(int start, int end) { 
		breakPoints = new ArrayList<>();
		breakPoints.add(new Range(start, end));
	}
	
	public BreakPoints() { }
	
	public void init(List<Integer> breakPoints) { 
		this.breakPoints = new ArrayList<>();
		for (int i = 0; i < breakPoints.size(); i=i+2)
			this.breakPoints.add(new Range(breakPoints.get(i), breakPoints.get(i+1)));
	}
	
	public BreakPoints(List<Range> breakPoints) { 
		this.breakPoints = new ArrayList<>(breakPoints);
	}

		
	
	/**
	 * computes new break points list based on a new break point introduced
	 * @param breakpoint
	 */
	public void computeLeftAndRight(int breakpoint) {
		List<Range> leftBreakPointsList = new ArrayList<>();
		List<Range> rightBreakPointsList = new ArrayList<>();
			
		
		int i=0;
		
		while (breakPoints.get(i).to <= breakpoint) {
			leftBreakPointsList.add(breakPoints.get(i));
			i++;
		}
		if (breakPoints.get(i).from <= breakpoint) {
			leftBreakPointsList.add(new Range(breakPoints.get(i).from, breakpoint));
			rightBreakPointsList.add(new Range(breakpoint+1, breakPoints.get(i).to));
			i++;
		}
		while (i < breakPoints.size()) {
			rightBreakPointsList.add(breakPoints.get(i));
			i++;
		}
		
		leftBreakPoints= new BreakPoints(leftBreakPointsList);
		rightBreakPoints= new BreakPoints(rightBreakPointsList);
	}

	/**
	 * gets all the info from left of a new breakpoint 
	 * @param position
	 * @return
	 */
	public BreakPoints getLeft() {
		return leftBreakPoints.copy();
	}
	
	/**
	 * gets all the info from the right of a new breakpoint 
	 * @param position
	 * @return
	 */
	public BreakPoints getRight() {
		return rightBreakPoints.copy();
	}

	/**
	 * gets the difference between the first and last event
	 * @return
	 */
	public double getLength() {
		return breakPoints.get(breakPoints.size()-1).to-breakPoints.get(0).from;
	}
	
	/**
	 * gets the total amount of genetic material encompassed
	 * by this break point as defined by 
	 * @return
	 */
	public double getGeneticLength() {
		int l = 0;
		for (int i = 0; i < breakPoints.size(); i++)
			l += breakPoints.get(i).size();
		return l;
	}
	
	public boolean isEmpty() {
		return breakPoints.size()==0;
	}
	

	public boolean withinLimits(int breakpoint) {
		if (this.breakPoints.get(0).from<=breakpoint && this.breakPoints.get(this.breakPoints.size()-1).to>breakpoint)
			return true;
		
		return false;
		
	} 
	
	public BreakPoints copy() {
		BreakPoints newBreakPoints = new BreakPoints(breakPoints);
		return newBreakPoints;
	}
		
	public String toString() {
		String val="";
		for (int i = 0; i < this.breakPoints.size(); i++) {
			val = val + "," + this.breakPoints.get(i).toString();
		}
				
		return val.substring(1);
		
	}
	
	
	
	public static class RangeComparator implements Comparator<Range> {
	    final int lessThan = -1;
	    final int greaterThan = 1;

	    @Override
	    public int compare(Range ra, Range rb) {
	    	if (ra.from<rb.from)
	    		return lessThan;
	    	else
	    		return greaterThan;
	    	
	    }

	}

	/**
	 * compute the intersection between this.breakpoonts and breakpoints
	 * @param breakPoints
	 */
	public void and(BreakPoints breakPoints) {
		List<Range> newBreaks = new ArrayList<>();
		int j = 0;
		

		for (int i = 0; i < this.breakPoints.size(); i++) {
			while (breakPoints.breakPoints.get(j).to <= this.breakPoints.get(i).from) {
				j++;
				if (j==breakPoints.breakPoints.size()) {
					this.breakPoints = new ArrayList<>(newBreaks);
					return;
				}

			}

			while (breakPoints.breakPoints.get(j).from <= this.breakPoints.get(i).to) {
				Range newR = this.breakPoints.get(i).getOverlap(breakPoints.breakPoints.get(j));
				if (newR!=null)
					newBreaks.add(newR);
				 
				if (breakPoints.breakPoints.get(j).to <= this.breakPoints.get(i).to)
					j++;
				else
					break;
				
				if (j==breakPoints.breakPoints.size()) {
					this.breakPoints = new ArrayList<>(newBreaks);
					return;
				}				
			}	
		}
		this.breakPoints = new ArrayList<>(newBreaks);
	}

	
	/**
	 * remove breakpoints from this.breakpoints
	 * @param breakPoints
	 */
	public void andNot(BreakPoints breakPoints) {
		List<Range> newBreaks = new ArrayList<>();
		int j = 0;
		

		
		for (int i = 0; i < this.breakPoints.size(); i++) {
			boolean rangeAdded = false;
			while (breakPoints.breakPoints.get(j).to <= this.breakPoints.get(i).from) {
				j++;
				if (j==breakPoints.breakPoints.size()) {
					if (!rangeAdded)
						newBreaks.add(this.breakPoints.get(i));
					i++;
					while (i < this.breakPoints.size()) {
						newBreaks.add(this.breakPoints.get(i));
						i++;
					}
					this.breakPoints = new ArrayList<>(newBreaks);
					return;
				}

			}

			while (breakPoints.breakPoints.get(j).from <= this.breakPoints.get(i).to) {
				List<Range> newR = this.breakPoints.get(i).getRemoved(breakPoints.breakPoints.get(j));
				if (newR!=null) {
					newBreaks.addAll(newR);
					rangeAdded = true;
				}
				 
				if (breakPoints.breakPoints.get(j).to <= this.breakPoints.get(i).to)
					j++;
				else
					break;
				
				if (j==breakPoints.breakPoints.size()) {
					if (!rangeAdded)
						newBreaks.add(this.breakPoints.get(i));
					i++;
					while (i < this.breakPoints.size()) {
						newBreaks.add(this.breakPoints.get(i));
						i++;
					}
					this.breakPoints = new ArrayList<>(newBreaks);
					return;
				}	
			}	
			
			if (!rangeAdded)
				newBreaks.add(this.breakPoints.get(i));
				
		}
		this.breakPoints = new ArrayList<>(newBreaks);
	}

	/**
	 * combine two lists of breakpoints
	 * @param breakPoints
	 */
	public void or(BreakPoints breakPoints) {
		// make a new list containing all breakpoints
		List<Range> newBreaks = new ArrayList<>();
		this.breakPoints.addAll(breakPoints.breakPoints);
		Collections.sort(this.breakPoints, rc);
				
		int nextto = this.breakPoints.get(0).to;
		int lastfrom = this.breakPoints.get(0).from;
		for (int i = 1; i < this.breakPoints.size(); i++) {
			if (this.breakPoints.get(i).from>(nextto+1)) {
				newBreaks.add(new Range(lastfrom, nextto));
				nextto = this.breakPoints.get(i).to;
				lastfrom = this.breakPoints.get(i).from;
			}else {
				nextto = Math.max(nextto, this.breakPoints.get(i).to);
			}
				
		}
		newBreaks.add(new Range(lastfrom, nextto));				
					
		this.breakPoints = new ArrayList<>(newBreaks);		
	}
	
	
	public Range getNewRange(int from, int to) {
		return new Range(from, to);
	}

	public class Range{
		public int from;
		public int to;
	
		public Range(int from, int to) {
			this.from = from;
			this.to = to;
		}
				
		public void combine(Range range) {
			this.from = Math.min(from, range.from);
			this.to = Math.max(to, range.to);
		}
		
		public int size() {
			return to-from+1;
		}
		
		public boolean isSmaller(Range range) {
			if (to < range.from)
				return true;
			return false;
		}
		
		public boolean isLarger(Range range) {
			if (from > range.to)
				return true;
			return false;
		}
		
		public Range getOverlap(Range range) {
			int newfrom = Math.max(from, range.from);
			int newto = Math.min(to, range.to);
			if (newto<newfrom)
				return null;
			
			return new Range(newfrom, newto);			
		}
		
		public List<Range> getRemoved(Range range) {
			// get the overlap between the two
			Range overlap = getOverlap(range);
			if (overlap==null)
				return null;
			
			List<Range> removed = new ArrayList<>();
			
			// remove the overlap from this.range
			if (from < overlap.from) {
				removed.add(new Range(from, overlap.from-1));
			}
			
			if (to > overlap.to) {
				removed.add(new Range(overlap.to+1, to));
			}
				
			return removed;
		}

		
		public String toString() {
			return from +"-"+to;
			
		}
	}

	
}
