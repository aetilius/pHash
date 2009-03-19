import java.util.*;

public class pHash
{
	native static BitSet videoHash(String file);
	native static int[]  audioHash(String file);
	native static BitSet imageHash(String file);
	static {
		System.loadLibrary("pHash-jni");
	}
}
