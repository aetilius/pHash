import java.util.*;

public class pHash
{
	native static long videoHash(String file);
	native static int[]  audioHash(String file);
	native static long imageHash(String file);
	native static int imageDistance(long hash1, long hash2);
	native static int audioDistance(int[] hash1, int[] hash2);
	static int videoDistance(long hash1, long hash2)
	{
		return imageDistance(hash1, hash2);
	}
	static {
		System.loadLibrary("pHash-jni");
	}


	public static void main(String args[])
	{
		long hash = imageHash(args[0]);
		long hash2 = imageHash(args[1]);
		System.out.println(hash);
		System.out.println(hash2);
		System.out.println(imageDistance(hash,hash2));

	}


}
