package org.phash.ni.image;

import static org.junit.Assert.*;
import org.junit.Test;

import java.io.File;

public class ImageHashTest {

	@Test
	public void testDCTImageHash(){
		ClassLoader classLoader = getClass().getClassLoader();
		File file = new File(classLoader.getResource("img1.jpg").getFile());
		assertTrue(file.exists());

		long hashval = ImageHash.dctImageHash(file.getPath());
		assertTrue(hashval != 0);
	}
	
}
