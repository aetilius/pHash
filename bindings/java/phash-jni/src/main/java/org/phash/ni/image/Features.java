package org.phash.ni.image;

public class Features {

	public double[] features;

	public Features(int length){
		this.features = new double[length];
	}

	public Features(double[] arr){
		this.features = arr;
	}

	public double[] getFeatures(){
		return features;
	}
}
