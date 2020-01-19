package org.phash.ni.image;

public class BMBHash {

	public byte[] hash;

	public BMBHash(int length){
		this.hash = new byte[length];
	}

	public BMBHash(byte[] arr){
		this.hash = arr;
	}

	public byte[] getHash(){
		return hash;
	}
}
