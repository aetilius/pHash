import java.util.*;

public class pHash
{
	class VideoMVPTree
	{
		native VideoMVPTree createTree(VideoHashRecord[] hashes);
		public VideoMVPTree(VideoHashRecord[] hashes)
		{
			createTree(hashes);
		}
	}
	native static long videoHash(String file);
	native static int[]  audioHash(String file);
	native static long imageHash(String file);
	native static int imageDistance(long hash1, long hash2);
	native static double audioDistance(int[] hash1, int[] hash2);
	static int videoDistance(long hash1, long hash2)
	{
		return imageDistance(hash1, hash2);
	}
	static {
		System.loadLibrary("pHash-jni");
	}


	public static void main(String args[])
	{
		//long hash = imageHash(args[0]);
		//long hash2 = imageHash(args[1]);
		//System.out.println(imageDistance(hash,hash2));
		int[] audioHash1 = audioHash(args[0]);
		
		int[] audioHash2 = audioHash(args[1]);
		System.out.println("Hashes match: " + Arrays.equals(audioHash1,audioHash2));
System.out.println("cs = " + audioDistance(audioHash1,audioHash2));
//		long vHash = videoHash(args[4]);
//		long vHash2 = videoHash(args[5]);
//		System.out.println(videoDistance(vHash,vHash2));

	}


}
