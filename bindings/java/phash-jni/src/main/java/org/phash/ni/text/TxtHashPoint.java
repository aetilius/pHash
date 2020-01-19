package org.phash.ni.text;

public class TxtHashPoint {

	public long hash;

	public long index;

	public TxtHashPoint(long hash, long index){
		this.hash = hash;
		this.index = index;
	}

	public long getHash(){
		return hash;
	}

	public long getIndex(){
		return index;
	}

	public void setDataMembers(long hash, long index){
		this.hash = hash;
		this.index = index;
	}
}
