using System;
using System.Runtime.InteropServices;

public class pHash
{
	static pHash() {
		Environment.SetEnvironmentVariable("PATH", Environment.GetEnvironmentVariable("PATH")+
							":/usr/local/lib");
	}
	[StructLayout(LayoutKind.Sequential)]
	public struct DP
	{
		public string id;
		public IntPtr hash;
		public IntPtr path;
		public UInt16 hash_length;
		public byte hash_type; 
	}
	public enum MVPRetCode
	{
		PH_SUCCESS = 0,   /* success */
    		PH_ERRPGSIZE,     /* page size error */
    		PH_ERRSMPGSIZE,
    		PH_ERRFILEOPEN,	  /* file open err */
    		PH_ERRFILETRUNC,
    		PH_ERRFILESEEK,
    		PH_ERRMMAP,        /* mmap'ing error */
    		PH_ERRMSYNC,       /* msync error */
    		PH_ERRTRUNC,       /* error truncating file */
    		PH_ERRSAVEMVP,      /* could not save mvp file */
    		PH_ERRARG,   /* general arg err*/
    		PH_ERRNULLARG,   /* null arg */
    		PH_ERRMEM,       /* general memory error  */
    		PH_ERRMEMALLOC,  /* mem alloc error */
    		PH_ERRNTYPE,      /* unrecognized node type */
    		PH_ERRCAP,     /* more results found than can be supported in ret array */
    		PH_ERRFILETYPE,  /*unrecognized file type  */
    		PH_ERRDIST,     
	}

	public enum HashType {
	    	BYTEARRAY   = 1,          /* refers to bitwidth of the hash value */
    		UINT16ARRAY = 2,
    		UINT32ARRAY = 4,
	 	UINT64ARRAY = 8,
	}

	[StructLayout(LayoutKind.Sequential)]
	public struct MVPFile
	{
		public string filename;
		private IntPtr buf;
		private int file_pos;
		private int fd;
		private byte filenumber;
		private byte nbdbfiles;
		public byte branchfactor;
		public byte pathlength;
		public byte leafcapacity;
		public int pgsize;
		public HashType hash_type;
		public HashCB hashdist;

		[UnmanagedFunctionPointer(CallingConvention.Cdecl)]
  		public delegate float HashCB(ref DP a, ref DP b);

	}
	[DllImport("pHash", CharSet=CharSet.Ansi)]
	private static extern MVPRetCode ph_query_mvptree(ref MVPFile m, ref DP query, int knearest,
						float radius, float threshold, IntPtr[] results,
						ref int count);

	public static MVPRetCode query_mvptree(ref MVPFile m, ref DP query, int knearest,
					float radius, float thresh, ref DP[] results, ref int count)
	{

		IntPtr[] points = new IntPtr[results.Length];
		MVPRetCode ret = ph_query_mvptree(ref m, ref query, knearest, radius, thresh, points,
							ref count);
		int i = 0;
		foreach(IntPtr ptr in points)
		{
			DP dp = (DP)Marshal.PtrToStructure(ptr, typeof(DP));
			results[i++] = dp;
			free(ptr);
		}
		return ret;
	}
	[DllImport("pHash", CharSet=CharSet.Ansi)]
	public static extern float ph_hamming_distance(ref DP a, ref DP b);

	[DllImport("pHash", CharSet=CharSet.Ansi)]
	private static extern IntPtr ph_malloc_datapoint(int t, int pl);

	[DllImport("pHash", CharSet=CharSet.Ansi)]
	public static extern void ph_mvp_init(ref MVPFile f);

	[DllImport("pHash", CharSet=CharSet.Ansi)]
	private static extern MVPRetCode ph_add_mvptree(ref MVPFile f, IntPtr[] points, int num, 
							ref int saved);
	public static MVPRetCode add_mvptree(ref MVPFile f, DP[] points, ref int numSaved)
	{
		IntPtr[] datapoints = new IntPtr[points.Length];
		int size = Marshal.SizeOf(typeof(DP));
		int i = 0;
		foreach(DP dp in points)
		{
			IntPtr ptr = Marshal.AllocHGlobal(size);
			Marshal.StructureToPtr(dp, ptr, false);
			datapoints[i++] = ptr;
		}
		MVPRetCode ret = ph_add_mvptree(ref f, datapoints, datapoints.Length, ref numSaved);
		for(i = 0; i < datapoints.Length; ++i)
			Marshal.FreeHGlobal(datapoints[i]);
		
		return ret;
		
	}
	
	[DllImport("pHash", CharSet=CharSet.Ansi)]
	public static extern int ph_dct_imagehash(string file, ref ulong hash);


	[DllImport("pHash", CharSet=CharSet.Ansi)]
	private static extern MVPRetCode ph_save_mvptree(ref MVPFile m, IntPtr[] points, int num);
	public static MVPRetCode save_mvptree(ref MVPFile m, DP[] points)
	{
		IntPtr[] datapoints = new IntPtr[points.Length];
		int size = Marshal.SizeOf(typeof(DP));
		int i = 0;
		foreach(DP dp in points)
		{
			IntPtr ptr = Marshal.AllocHGlobal(size);
			Marshal.StructureToPtr(dp, ptr, false);
			datapoints[i++] = ptr;
			
		}
		MVPRetCode ret = ph_save_mvptree(ref m, datapoints, datapoints.Length);
		for(i = 0; i < datapoints.Length; ++i)
			Marshal.FreeHGlobal(datapoints[i]);

		return ret;
	}
	
	[DllImport("pHash", CharSet=CharSet.Ansi)]
	public static extern IntPtr ph_mh_imagehash(string file, ref int n, float alpha,
						float lvl);

	[DllImport("pHash", CharSet=CharSet.Ansi)]
	public static extern IntPtr ph_audiohash(IntPtr buf, int nbbuf, UInt32[] hashbuf,
						int nbcap, int sr, ref int nbframes);

	[DllImport("pHash", CharSet=CharSet.Ansi)]
	private static extern IntPtr ph_readaudio(string file, int sr, int ch, 
					float[] buf, ref int buflen,
					float nbsecs);
	[DllImport("libc")]
	public static extern void free(IntPtr p);
	
}
