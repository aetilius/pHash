#include "pHash.h"

static int _bmb_setbit(BMBHash &bh, uint32_t bit) {
	if (bh.hash == NULL) return -1;
	bh.hash[bit/8] |= 0x01 << (bit%8);
	return 0;
}
static void _ph_bmb_new(BMBHash &bh, uint32_t n_bits) {
	bh.bytelength = n_bits/8 + 1;
	bh.hash = new uint8_t[bh.bytelength];
	for (int i=0;i<bh.bytelength;i++){
		bh.hash[i] = 0x00;
	}
}

void ph_bmb_free(BMBHash &bh) {
	delete[] bh.hash;
}

int ph_bmb_imagehash(const char *file, BMBHash &ret_hash) {
	if (!file) return -1;

	CImg<uint8_t> src;
	const float sigma = 1.0f;
	const int preset_width = 256;
	const int preset_height = 256;
	const int blk_width = 16;
	const int blk_height = 16;
	int n_blocks = (preset_width*preset_height)/(blk_width*blk_height); 
	
	try {
		src.load(file);
	} catch (CImgIOException ex) {
		return -1;
	}

	CImg<float> img; 
	switch (src.spectrum()) {
	case 4:
	case 3: // from RGB
		img = src.get_RGBtoYUV();
		img.channel(0);
		break;
	case 1: // grayscale
		break;;
	default:
		return -1;
	}
	
	img.blur(sigma, false, true);  /* gaussian blur with dirichlet boundary condition */
	img.resize(preset_width, preset_height, 1, 1, 1); /* linear interpolation */

	CImg<double> localmeans(n_blocks, 1, 1, 1, 0);
	
	int blockidx = 0;
	for (int blockrow = 0; blockrow < preset_height - blk_height; blockrow += blk_height) {
		for (int blockcol = 0; blockcol < preset_width - blk_width; blockcol += blk_width) {
			float acc = 0;
			for (int subrow=blockrow;subrow < blockrow+blk_height;subrow++){
				for (int subcol=blockcol;subcol < blockcol+blk_width;subcol++){
					float pxl = img.atXY(subcol, subrow);
					acc = acc + pxl;
				}
			}
			localmeans(blockidx) = acc / (blk_height*blk_height);
			blockidx++;
		}
	}

	double median_value = localmeans.median();
	
	_ph_bmb_new(ret_hash, n_blocks);

	for (int i = 0; i < n_blocks; i++) {
		if (localmeans(i) > median_value) {
			_bmb_setbit(ret_hash, i);
		} 
	}
   	return 0;
}

int ph_bmb_distance(const BMBHash &bh1, const BMBHash &bh2){

	return ph_hammingdistance2(bh1.hash, bh1.bytelength, bh2.hash, bh2.bytelength);
}

