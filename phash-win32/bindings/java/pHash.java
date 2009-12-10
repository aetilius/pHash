	abstract class Hash
	{
		String filename;
	}
	class AudioHash extends Hash 
	{
		int[] hash;
	}
	class TextHash extends Hash 
	{
		String[] hash;
	}
	class VideoHash extends Hash 
	{
		long hash;
	}
	class ImageHash extends Hash 
	{
		long hash;
	}

class pHash
{
	class MVPTree
	{
		String mvpFile;

		public MVPTree(String filename) { mvpFile = filename; }
		public native boolean create(Hash[] hashes);
		public native Hash[] query(Hash hash, float radius, int maxResults);
		public native boolean add(Hash[] hashes);
	}

	native static VideoHash videoHash(String file);
	native static AudioHash audioHash(String file);
	native static ImageHash imageHash(String file);
	native static TextHash textHash(String file);
	native static int imageDistance(ImageHash hash1, ImageHash hash2);
	native static double audioDistance(AudioHash hash1, AudioHash hash2);
	native static int videoDistance(VideoHash hash1, VideoHash hash2);
	native static int textDistance(TextHash txtHash1, TextHash txtHash2);
	private native static void pHashInit();
	private native static void cleanup();
	static {
		System.loadLibrary("pHash-jni");
		pHashInit();
	}


	public static void main(String args[])
	{
			int i = 0;
			if(args[i].equals("-a"))
			{
				AudioHash audioHash1 = audioHash(args[1]);
				AudioHash audioHash2 = audioHash(args[2]);
				System.out.println("cs = " + audioDistance(audioHash1,audioHash2));
			}
			else if(args[i].equals("-i"))
			{
				ImageHash imHash = imageHash(args[1]);
				ImageHash imHash2 = imageHash(args[2]);
				System.out.println("File 1: " + imHash.filename);
				System.out.println("File 2: " + imHash2.filename);

				System.out.println(imageDistance(imHash,imHash2));
			}
			else if(args[i].equals("-v"))
			{
				VideoHash vHash = videoHash(args[1]);
				VideoHash vHash2 = videoHash(args[2]);
				System.out.println(videoDistance(vHash,vHash2));
			}
			else if(args[i].equals("-t"))
			{
				TextHash txtHash = textHash(args[1]);
				TextHash txtHash2 = textHash(args[2]);
                                System.out.println(textDistance(txtHash,txtHash2));
			}

			pHash.cleanup();

	}


}
