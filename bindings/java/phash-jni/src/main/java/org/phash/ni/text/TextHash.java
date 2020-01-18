package org.phash.ni.text;

public class TextHash {
	static public native TxtHashPoint[] textHash(String file);
	
	static public native TxtMatch[] compareTextHashes(TxtHashPoint[] hash1, TxtHashPoint[] hash2);

	static {
		System.loadLibrary("pHash-jni");
	}
}
