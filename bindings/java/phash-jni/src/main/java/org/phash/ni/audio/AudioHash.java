package org.phash.ni.audio;

public class AudioHash {
	static public native int countAudioSamples(String file, int sr, int nChannels);
	static public native float[] readAudio(String file, int sr, int nChannels, float nbSeconds);
	static public native int[] audioHash(float[] buf, int sr);
	static public native double[] audioDistanceBER(int[] hash1, int[] hash2, float threshold, int blocksize);

	static {
		System.loadLibrary("pHash-jni");
	}
}
