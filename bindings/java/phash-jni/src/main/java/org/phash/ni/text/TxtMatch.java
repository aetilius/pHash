package org.phash.ni.text;

public class TxtMatch {

	public long firstIndex;

	public long secondIndex;

	public long length;

	public TxtMatch(long firstIndex, long secondIndex, long length){
		this.firstIndex = firstIndex;
		this.secondIndex = secondIndex;
		this.length = length;
	}

	public long getFirstIndex(){
		return firstIndex;
	}

	public long getSecondIndex(){
		return secondIndex;
	}

	public long getLengthIndex(){
		return length;
	}

	public void setDataMembers(long firstIndex, long secondIndex, long length){
		this.firstIndex = firstIndex;
		this.secondIndex = secondIndex;
		this.length = length;
	}
	
}
