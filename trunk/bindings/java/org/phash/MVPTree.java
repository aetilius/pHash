package org.phash;
public class MVPTree
{
	protected String mvpFile;
	public MVPTree(String filename) { 
		mvpFile = filename; 
	}
	public native boolean create(Hash[] hashes);
	public native Hash[] query(Hash hash, float radius, float thresh,
				int maxResults);
	public native boolean add(Hash[] hashes);
}

