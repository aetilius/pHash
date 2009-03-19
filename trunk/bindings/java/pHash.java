import java.util.*;

public class pHash
{
	native static BitSet calculateHash(String file, int hashType);
	static {
		System.loadLibrary("pHash-jni");
	}
}
