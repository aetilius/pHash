package org.phash.ni.image;

public class Digest {

	public byte[] coeffs;

	public Digest(int length){
		this.coeffs = new byte[length];
	}

	public Digest(byte[] arr){
		this.coeffs = arr;
	}

	public byte[] getDigest(){
		return coeffs;
	}
}
