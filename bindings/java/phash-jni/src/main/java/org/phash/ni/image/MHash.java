package org.phash.ni.image;

public class MHash {

	public byte[] hash;

	public MHash(int length){
		this.hash = new byte[length];
	}

	public MHash(byte[] arr){
		this.hash = arr;
	}

	public byte[] getHash(){
		return hash;
	}
}
