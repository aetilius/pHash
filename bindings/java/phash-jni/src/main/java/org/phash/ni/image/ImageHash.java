package org.phash.ni.image;

public class ImageHash {

	static public native long dctImageHash(String file);

	static public native int dctImageDistance(long hash1, long hash2);
	
	static public native BMBHash bmbImageHash(String file);

	static public native int bmbDistance(BMBHash hash1, BMBHash hash2);

	static public native Features radialFeatureVector(String file, int nAngles, double sigma, double gamma);
	
	static public native Digest imageDigest(String file, double sigma, double gamma, int nAngles);

	static public native double peakCrossCorr(Digest digest1, Digest digest2);
	
	static public native MHash mhImageHash(String file, float alpha, float levl);

	static public native int mhImageDistance(MHash hash1, MHash hash2);

	static {
		System.loadLibrary("pHash-jni");
	}

	
}
