import java.io.*;

public class pHash
{

	native static VideoHash videoHash(String file);
	native static AudioHash audioHash(String file);
	native static DCTImageHash dctImageHash(String file);
	native static MHImageHash mhImageHash(String file);
	native static TextHash textHash(String file);
	native static double imageDistance(ImageHash hash1, ImageHash hash2);
	native static double audioDistance(AudioHash hash1, AudioHash hash2);
	native static double videoDistance(VideoHash hash1, VideoHash hash2, int threshold);
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
			if(args[i].equals("-mvp"))
			{
				File dir = new File(args[1]);
				if(dir.isDirectory())
				{
					File[] files = dir.listFiles();
					MHImageHash[] hashes = new MHImageHash[files.length];
					for(int j = 0; j < files.length; ++j)
					{
						hashes[j] = mhImageHash(files[j].toString());
												
					}	
					MVPTree mvp = new MVPTree("mvp");
					boolean result = mvp.create(hashes);
					if(result)
						System.out.println("Successfully created MVP tree");									
				}
			}
			else if(args[i].equals("-a"))
			{
				AudioHash audioHash1 = audioHash(args[1]);
				AudioHash audioHash2 = audioHash(args[2]);
				System.out.println("cs = " + audioDistance(audioHash1,audioHash2));
			}
			else if(args[i].equals("-dct"))
			{
				DCTImageHash imHash = dctImageHash(args[1]);
				DCTImageHash imHash2 = dctImageHash(args[2]);
				System.out.println("File 1: " + imHash.filename);
				System.out.println("File 2: " + imHash2.filename);

				System.out.println(imageDistance(imHash,imHash2));
			}
			else if(args[i].equals("-mh"))
			{
				MHImageHash imHash = mhImageHash(args[1]);
				MHImageHash imHash2 = mhImageHash(args[2]);
				System.out.println("File 1: " + imHash.filename);
				System.out.println("File 2: " + imHash2.filename);

				System.out.println(imageDistance(imHash,imHash2));

			}
			else if(args[i].equals("-v"))
			{
				VideoHash vHash = videoHash(args[1]);
				VideoHash vHash2 = videoHash(args[2]);
				System.out.println(videoDistance(vHash,vHash2, 21));
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
