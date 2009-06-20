#include <stdio.h>
#include <dirent.h>
#include "pHash.h"


DP* ph_read_datapoint(MVPFile *m){
    DP *dp = NULL;
    uint8_t active;
    uint16_t byte_len;
    uint16_t id_len;
    uint16_t hash_len;
    int type = m->hash_type;
    int PathLength = m->pathlength;
    off_t file_pos = m->file_pos;
    off_t offset_mask = m->internal_pgsize - 1; /*BUG: how do i know which page size ??? */

    memcpy(&active, &(m->buf[file_pos & offset_mask]), sizeof(uint8_t));
    file_pos++;

    memcpy(&byte_len,&(m->buf[file_pos & offset_mask]), sizeof(uint16_t));
    file_pos += sizeof(uint16_t);

    if ((active == 0) &&(byte_len == 0)){
	return dp;
    }
    dp = ph_malloc_datapoint(type,PathLength);

    memcpy(&id_len,&(m->buf[file_pos & offset_mask]), sizeof(uint16_t));
    file_pos += sizeof(uint16_t);

    dp->id = (char*)malloc((id_len+1)*sizeof(uint8_t));
    memcpy(dp->id, &(m->buf[file_pos & offset_mask]), id_len);
    dp->id[id_len] = '\0';
    file_pos += id_len;

    memcpy(&hash_len, &(m->buf[file_pos & offset_mask]), sizeof(uint16_t));
    dp->hash_length = hash_len;
    file_pos += sizeof(uint16_t);

    dp->hash = malloc(hash_len*type);
    memcpy(dp->hash, &(m->buf[file_pos & offset_mask]), hash_len*type);
    file_pos += hash_len*type;
    
    memcpy(dp->path, &(m->buf[file_pos & offset_mask]), PathLength*sizeof(float));
    file_pos += PathLength*sizeof(float);

    m->file_pos = file_pos;

    return dp;
}


MVPFile* _ph_map_mvpfile(uint8_t filenumber, off_t offset, MVPFile *m){
    MVPFile *ret_file = NULL;

    if (filenumber == 0){ /* in the file denoted by m , must advance to page containing offset */
	off_t page_mask = ~(m->internal_pgsize - 1);
	off_t page_offset = offset & page_mask;
	ret_file = m;
	m->isleaf = 0;
	m->file_pos = offset;
	if (munmap(m->buf, m->internal_pgsize) < 0){
	    perror("munmap");
	}
	m->buf = (char*)mmap(NULL,m->internal_pgsize,PROT_WRITE|PROT_READ,MAP_SHARED,m->fd, page_offset);
	if (m->buf == MAP_FAILED){
	    perror("mmap");
	    ret_file = NULL;
	}
	
    } else { /* open and map to new file denoted by m->filename and filenumber */
	off_t page_mask = ~(m->leaf_pgsize - 1);
	ret_file = (MVPFile*)malloc(sizeof(MVPFile));
        if (!ret_file)
	    return NULL;
        char extfile[256];
	sprintf(extfile, "%s%d.mvp", m->filename, filenumber);
        ret_file->filename = strdup(m->filename);
	ret_file->fd = open(extfile, O_RDWR);
	ret_file->pathlength = m->pathlength;
	ret_file->hash_type = m->hash_type;
	ret_file->hashdist = m->hashdist;
	ret_file->branchfactor = m->branchfactor;
	ret_file->leafcapacity = m->leafcapacity;
	ret_file->nbdbfiles  = m->nbdbfiles;
	ret_file->internal_pgsize = m->internal_pgsize;
	ret_file->leaf_pgsize = m->leaf_pgsize;
	ret_file->isleaf = 1;
	ret_file->file_pos = offset;
	off_t page_offset = offset & page_mask;
	
	ret_file->buf = (char*)mmap(NULL, ret_file->leaf_pgsize, PROT_READ|PROT_WRITE, MAP_SHARED, 
                                         ret_file->fd, page_offset);
	if (ret_file->buf == MAP_FAILED){
	    perror("mmap");
	    return NULL;
	}

    }
    return ret_file;
}

void _ph_unmap_mvpfile(uint8_t filenumber, off_t orig_pos, MVPFile *m, MVPFile *m2){

    if (filenumber == 0){ /*remap to same main file  */
	off_t page_mask = ~(m->internal_pgsize - 1);
	if (munmap(m->buf, m->internal_pgsize) < 0){
	    perror("munmap");
	}
	m->file_pos = orig_pos;
	off_t pg_offset = orig_pos & page_mask;
	m->buf = (char*)mmap(NULL, m->internal_pgsize, PROT_WRITE|PROT_READ, MAP_SHARED,m->fd,pg_offset);
	if (m->buf == MAP_FAILED){
	    perror("mmap");
	}
    } else { /*remap to back to main file  */
	if (munmap(m2->buf, m2->leaf_pgsize) < 0){
	    perror("munmap");
	}
	if (close(m2->fd)<0)
	    perror("close");
        free(m2);
    }
}


MVPRetCode ph_query_mvptree(MVPFile *m, DP *query, int knearest, float radius, 
                                              DP **results, int *count, int level){
    int BranchFactor = m->branchfactor;
    int LengthM1 = BranchFactor-1;
    int LengthM2 = BranchFactor*LengthM1;
    int PathLength = m->pathlength;
    int res = PH_SUCCESS;
    if ((!m)||(!query))
	return PH_ERR_ARGLIST;

    hash_compareCB hashdist = m->hashdist;

    if (!hashdist)
	return PH_ERR_NODISTFUNC;

    off_t offset_mask, page_mask;
    if (m->isleaf){
	offset_mask = m->leaf_pgsize - 1;
	page_mask = ~(m->leaf_pgsize - 1);
    }
    else {
	offset_mask = m->internal_pgsize - 1;
	page_mask = ~(m->internal_pgsize - 1);
    }

    uint8_t ntype;
    memcpy(&ntype, &m->buf[m->file_pos & offset_mask], sizeof(uint8_t));
    m->file_pos++;

    if (ntype == 0){ /* leaf */
	DP *sv1 = ph_read_datapoint(m);
	DP *sv2 = ph_read_datapoint(m);
	float d1 = hashdist(query,sv1);

	/* check if distance(sv1,query) <= radius  */
	if (d1 <= radius){
	    results[(*count)++] = sv1;
	    if (*count >= knearest)
		return PH_RESULTSFULL;
	} else {
	    ph_free_datapoint(sv1);
	}

	if (sv2){
	    float d2 = hashdist(query,sv2);
	    /* check if distance(sv2,query) <= radius */
	    if (d2 <= radius){
		results[(*count)++] = sv2;
		if (*count >= knearest)
		    return PH_RESULTSFULL;
	    } else {
		ph_free_datapoint(sv2);
	    }

	    uint8_t Np;
	    memcpy(&Np, &m->buf[m->file_pos & offset_mask], sizeof(uint8_t));
	    m->file_pos += sizeof(uint8_t);
	    off_t curr_pos;

	    /* read in each datapoint in the leaf - only retrieve the point if it 
	       dist(sv1,dp)=da  and dist(sv2,dp)=db cannot preclude the point */
	    for (int i=0;i<Np;i++){
		int include = 1;
		float da, db;
		off_t point_offset;
		memcpy(&da, &m->buf[m->file_pos & offset_mask], sizeof(float));
		m->file_pos += sizeof(float);
		memcpy(&db, &m->buf[m->file_pos & offset_mask], sizeof(float));
		m->file_pos += sizeof(float);
		if ((d1-radius <= da)&&(d1+radius >= da)&&(d2-radius <= db)&&(d2+radius >= db)){
		    memcpy(&point_offset, &m->buf[m->file_pos & offset_mask],sizeof(off_t));
		    m->file_pos += sizeof(off_t);
		    curr_pos = m->file_pos;
		    m->file_pos = point_offset;
		    DP *dp = ph_read_datapoint(m);
		    m->file_pos = curr_pos;
		    
		    /* test each path[] distance and as soon as one does not fit 
		       disclude the point                                        */
                    if (dp){
			for (int j=0;j<level;j++){
			    if ((query->path[j]-radius <= dp->path[j])
				&&(query->path[j]+radius >= dp->path[j])){
				if (hashdist(query,dp) > radius){
				    include = 0;
				    break;
				}
			    }
			}
		    } 
		    if (include){
			results[(*count)++] = dp;
			if (*count >= knearest)
			    return PH_RESULTSFULL;
		    } else {
			ph_free_datapoint(dp);
		    }
		} else {
		    m->file_pos += sizeof(off_t);
		}
	    }
	}
    } else if (ntype == 1) { /* internal */
	/* read sv1, sv2 */
	DP *sv1 = ph_read_datapoint(m);
	DP *sv2 = ph_read_datapoint(m);
	/* read 1st and 2nd level pivots */
	float *M1 = (float*)malloc(LengthM1*sizeof(float));
	if (!M1){
	    return PH_MEMALLOC;
	}
	float *M2 = (float*)malloc(LengthM2*sizeof(float));
	if (!M2){
	    return PH_MEMALLOC;
	}

	memcpy(M1, &m->buf[m->file_pos & offset_mask], LengthM1*sizeof(float));
        m->file_pos += LengthM1*sizeof(float);
	memcpy(M2, &m->buf[m->file_pos & offset_mask], LengthM2*sizeof(float));
	m->file_pos += LengthM2*sizeof(float);

	float d1 = hashdist(query, sv1);
	float d2 = hashdist(query, sv2);

	/* fill in path values in query */
	if (level < PathLength)
	    query->path[level] = d1;
	if (level < PathLength - 1)
	    query->path[level+1] = d2;

	/* check if sv1 sv2 are close enough to query  */
	if (d1 <= radius){
	    results[(*count)++] = sv1;
	    if (*count >= knearest)
		return PH_RESULTSFULL;
	} else {
	    ph_free_datapoint(sv1);
	}

	if (d2 <= radius){
	    results[(*count)++] = sv2;
	    if (*count >= knearest)
		return PH_RESULTSFULL;
	} else {
	    ph_free_datapoint(sv2);
	}

	/* based on d1,d2 values, find appropriate child nodes to explore */
	
	int pivot1, pivot2;
	uint8_t filenumber;
	off_t child_pos, curr_pos;
	off_t start_pos = m->file_pos;
	off_t orig_pos;
	/* check <= each M1 pivot */
	for (pivot1=0;pivot1 < LengthM1;pivot1++){
	    if (d1-radius <= M1[pivot1]){
		/* check <= each M2 pivot */
		for (pivot2=0;pivot2<LengthM1;pivot2++){
		    if (d2 - radius <= M2[pivot2+pivot1*LengthM1]){
			/*determine pos from which to read filenumber and offset */
			curr_pos = start_pos + (pivot2+pivot1*BranchFactor)*(sizeof(uint8_t)+sizeof(off_t));
			m->file_pos = curr_pos;
			memcpy(&filenumber,&(m->buf[m->file_pos&offset_mask]),sizeof(uint8_t));
                        m->file_pos++;
			memcpy(&child_pos,&(m->buf[m->file_pos&offset_mask]),sizeof(off_t));
                        m->file_pos += sizeof(off_t);

			/*save position and remap to new file/position  */
			orig_pos = m->file_pos;
			MVPFile *m2 = _ph_map_mvpfile(filenumber,child_pos, m);
			if (m2){
			   res=ph_query_mvptree(m2,query,knearest,radius, results, count, level+2);
			}

			/* unmap and remap to the origional file/posion */
			_ph_unmap_mvpfile(filenumber, orig_pos, m, m2);
			
		    }
		}
		/* check > last M2 */
		if (d2+radius >= M2[LengthM1-1+pivot1*LengthM1]){

		    /*determine position from which to read filenumber and offset */
		    curr_pos = start_pos + (BranchFactor-1+pivot1*BranchFactor)*(sizeof(uint8_t)+sizeof(off_t));
		    m->file_pos = curr_pos;
		    memcpy(&filenumber,&m->buf[m->file_pos&offset_mask],sizeof(uint8_t));
		    m->file_pos++;
		    memcpy(&child_pos,&m->buf[m->file_pos&offset_mask],sizeof(off_t));
		    m->file_pos += sizeof(off_t);

		    /*saveposition and remap to new file/position */
                    orig_pos = m->file_pos;

		    MVPFile *m2 = _ph_map_mvpfile(filenumber, child_pos,m); 
		    if (m2){
			res = ph_query_mvptree(m2,query,knearest,radius,results,count,level+2);
		    }
		    /*unmap and remap to original file/position  */
		    _ph_unmap_mvpfile(filenumber, orig_pos, m, m2);
		}
	    }
	}
	/* check >=  last M1 pivot */
	if (d1+radius >= M1[LengthM1-1]){

	    /* check <= each M2 pivot */
	    for (pivot2=0;pivot2<LengthM1;pivot2++){
		if (d2-radius <= M2[pivot2+LengthM1*LengthM1]){

		    /*determine pos from which to read filenumber and position  */
		    curr_pos = start_pos + (pivot2+LengthM1*BranchFactor)*(sizeof(uint8_t)+sizeof(off_t));
		    m->file_pos = curr_pos;
		    memcpy(&filenumber,&m->buf[m->file_pos&offset_mask],sizeof(uint8_t));
		    m->file_pos++;
		    memcpy(&child_pos,&m->buf[m->file_pos&offset_mask],sizeof(off_t));
		    m->file_pos += sizeof(off_t);

		    /*save file position and remap to new filenumber/offset  */
		    orig_pos = m->file_pos;
		    MVPFile *m2 = _ph_map_mvpfile(filenumber, child_pos, m);
		    if (m2){
			res =ph_query_mvptree(m2, query, knearest, radius,results, count, level+2);
		    }
		    /* unmap/remap to original filenumber/position */
		    _ph_unmap_mvpfile(filenumber, orig_pos, m, m2);
		}
	    }

	    /* check >= last M2 pivot */
	    if (d2+radius >= M2[LengthM1-1+LengthM1*LengthM1]){

		/* determine position from which to read filenumber and child position */
		curr_pos = start_pos + (BranchFactor-1+LengthM1*BranchFactor)*(sizeof(uint8_t)+sizeof(off_t));
		m->file_pos = curr_pos;
		memcpy(&filenumber,&m->buf[m->file_pos&offset_mask],sizeof(uint8_t));
		m->file_pos += sizeof(uint8_t);
		memcpy(&child_pos,&m->buf[m->file_pos&offset_mask],sizeof(off_t));
		m->file_pos += sizeof(off_t);

		/* save position and remap to new filenumber/position */
		orig_pos = m->file_pos;
		MVPFile *m2 = _ph_map_mvpfile(filenumber, child_pos, m);
		if (m2){
		    res = ph_query_mvptree(m2, query, knearest, radius, results, count, level+2);
		}
		/* return to original and remap to original filenumber/position */
		_ph_unmap_mvpfile(filenumber, orig_pos, m, m2);
	    }
	}
} else { /* unrecognized node */
	return PH_ERR_NTYPE;
    }
    return (MVPRetCode)res;
}


MVPRetCode ph_query_mvptree(MVPFile *m, DP *query, int knearest, float radius, 
                                               DP **results, int *count){

    if (m->internal_pgsize == 0)
	m->internal_pgsize = sysconf(_SC_PAGE_SIZE);
    if (m->leaf_pgsize == 0)
	m->leaf_pgsize = sysconf(_SC_PAGE_SIZE);

    
    if ((m->internal_pgsize < sysconf(_SC_PAGE_SIZE))||(m->leaf_pgsize < sysconf(_SC_PAGE_SIZE)))
	return PH_ERRPGSIZE;

    /* pg size must be power of 2 */
    if (m->internal_pgsize & (m->internal_pgsize - 1))
	return PH_ERRPGSIZE;
    if (m->leaf_pgsize & (m->leaf_pgsize - 1))
	return PH_ERRPGSIZE;

    char mainfile[256];
    sprintf(mainfile, "%s.mvp", m->filename);
    m->fd = open(mainfile, O_RDWR);
    if (m->fd < 0){
	perror("fopen");
	return PH_ERRFILE;
    }
    m->file_pos = 0;
    m->buf=(char*)mmap(NULL,m->internal_pgsize,PROT_READ|PROT_WRITE,MAP_SHARED,m->fd,m->file_pos);
    if (m->buf == MAP_FAILED){
	perror("mmap");
	return PH_ERRMAP;
    }


    char tag[17];
    int version;
    int fileintpgsize, fileleafpgsize;
    uint8_t nbdbfiles, bf, p, k, type;

    memcpy((char*)tag, (char*)&(m->buf[m->file_pos]), 16);
    tag[16] = '\0';
    m->file_pos += 16;

    memcpy(&version, &m->buf[m->file_pos], sizeof(int));
    m->file_pos += sizeof(int);

    memcpy(&fileintpgsize, &m->buf[m->file_pos], sizeof(int));
    m->file_pos += sizeof(int);

    memcpy(&fileleafpgsize, &m->buf[m->file_pos], sizeof(int));
    m->file_pos += sizeof(int);

    memcpy(&nbdbfiles, &m->buf[m->file_pos++], 1);

    memcpy(&bf, &m->buf[m->file_pos++], 1);

    memcpy(&p, &m->buf[m->file_pos++], 1);

    memcpy(&k, &m->buf[m->file_pos++], 1);
    
    memcpy(&type, &m->buf[m->file_pos++], 1);

    if ((fileintpgsize != m->internal_pgsize) || (fileleafpgsize != m->leaf_pgsize))
	return PH_ERRPGSIZE;


    m->file_pos = HeaderSize;

    /* finish the query by calling the recursive auxiliary function */
    *count = 0;
    MVPRetCode res = ph_query_mvptree(m,query,knearest,radius,results,count,0);

    if (munmap(m->buf, m->internal_pgsize) < 0)
	perror("munmap");
    m->buf = NULL;

    if (close(m->fd) < 0)
	perror("fclose");
    m->fd = 0;
    m->file_pos = 0;

    return res;
}

float hammingdistance(DP *pntA, DP *pntB){

    uint8_t htypeA = pntA->hash_type;
    uint8_t htypeB = pntB->hash_type;
    if (htypeA != htypeB)
	return -1.0;
    if (htypeA != UINT64ARRAY)
	return -1.0;
    if ((pntA->hash_length > 1) || (pntB->hash_length > 1))
	return -1.0;
    ulong64 *hashA = (ulong64*)pntA->hash;
    ulong64 *hashB = (ulong64*)pntB->hash;
    int res = ph_hamming_distance(*hashA, *hashB);
    return (float) res;
}



int main(int argc, char **argv){
 
    const char *filename = argv[1];/* name of file to save db */

    MVPFile mvpfile;
    ph_mvp_init(&mvpfile);
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = hammingdistance;
    mvpfile.hash_type = UINT64ARRAY;


    ulong64 tmphash = 0;

    char c = 'y';
    char queryfile[256];
    float radius = 21.0;
    const int knearest = 20;
    DP *results[knearest];
    int nbfound;

    DP *query = ph_malloc_datapoint(mvpfile.hash_type, mvpfile.pathlength);
    do {
	nbfound = 0;
	printf("enter file:");
	fgets(queryfile, 100, stdin);
	queryfile[strlen(queryfile)-1] = '\0';
	if (ph_dct_imagehash(queryfile, tmphash) < 0)
	    continue;
	query->id = strdup(queryfile);
	query->hash = malloc(UINT64ARRAY);
	ulong64 *ptrHash = (ulong64*)query->hash;
	ptrHash[0] = tmphash; 

	printf("search for: %s\n", queryfile);
	int res;
	if ((res = ph_query_mvptree(&mvpfile, query, knearest, radius, results, &nbfound)) != 0){
	    printf("could not complete query, ret result %d\n",res);
	    continue;
	}
	printf(" %d files found\n", nbfound);
	for (int i=0;i<nbfound;i++){
	    ulong64 *hash = (ulong64*)results[i]->hash;
	    printf(" %d %s %llx\n", i, results[i]->id, *hash);
	}
	printf("again? y/n:");
	c = fgetc(stdin);
	getchar();
    } while (c != 'n');


    return 0;
}
