import java.util.*;

public class pHash
{
	native static long videoHash(String file);
	native static int[]  audioHash(String file);
	native static long imageHash(String file);
	native static int distance(long hash1, long hash2);
	static {
		System.loadLibrary("pHash-jni");
	}
}
