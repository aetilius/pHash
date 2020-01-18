package org.phash.ni.video;

public class VideoHash {

	static public native long[] dctVideoHash(String file);

	static public native double dctVideoHashDistance(long[] hash1, long[] hash2, int threshold);

	static {
		System.loadLibrary("pHash-jni");
	}
}
