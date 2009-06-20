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

    if ((active == 0) ||(byte_len == 0)){
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


off_t ph_save_datapoint(DP *dp, MVPFile *m){
    uint8_t active = 1;
    uint16_t byte_len = 0;
    off_t point_pos = m->file_pos;
    off_t offset_mask = m->internal_pgsize - 1; /*bug: how do i know which page size?   */
    if (dp == NULL){
	active = 0;
        memcpy(&(m->buf[m->file_pos & offset_mask]),&active, 1);
	m->file_pos++;
	memcpy(&(m->buf[m->file_pos & offset_mask]),&byte_len,sizeof(uint16_t));
	m->file_pos += sizeof(uint16_t);
	return point_pos;
    }
    int type = m->hash_type;
    int PathLength = m->pathlength;
    uint16_t id_len = strlen(dp->id);
    uint16_t hash_len = dp->hash_length;
    byte_len = id_len + 4 + hash_len*(m->hash_type) + PathLength*sizeof(float);

    memcpy(&(m->buf[m->file_pos & offset_mask]), &active, 1);
    m->file_pos++;

    memcpy(&(m->buf[m->file_pos & offset_mask]), &byte_len, sizeof(uint16_t));
    m->file_pos += sizeof(uint16_t);

    memcpy(&(m->buf[m->file_pos & offset_mask]), &id_len, sizeof(uint16_t));
    m->file_pos += sizeof(uint16_t);

    memcpy(&(m->buf[m->file_pos & offset_mask]), (dp->id), id_len);
    m->file_pos += id_len;

    memcpy(&(m->buf[m->file_pos & offset_mask]), &hash_len, sizeof(uint16_t));
    m->file_pos += sizeof(uint16_t);

    memcpy(&(m->buf[m->file_pos & offset_mask]), (dp->hash), hash_len*type);
    m->file_pos += hash_len*type;

    memcpy(&(m->buf[m->file_pos & offset_mask]), (dp->path), sizeof(float)*PathLength);

    m->file_pos += PathLength*sizeof(float);


    return point_pos;
}


FileIndex* ph_save_mvptree(MVPFile *m, DP **points, int nbpoints, int saveall_flag, int level){
    int Np = (nbpoints >= 2) ? nbpoints - 2 : 0; 
    int BranchFactor = m->branchfactor;
    int PathLength = m->pathlength;
    int LeafCapacity = m->leafcapacity;
    int LengthM1 = BranchFactor-1;
    int LengthM2 = BranchFactor*LengthM1;
    int Fanout = BranchFactor*BranchFactor;

    if ((!m) || (!points) || (nbpoints <= 0))
	return NULL;
    FileIndex *pOffset = (FileIndex*)malloc(sizeof(FileIndex));
    if (!pOffset){
        fprintf(stderr, "unable to allocate mem\n");
	return NULL;
    }
    hash_compareCB hashdist = m->hashdist;
    if (!hashdist){
	fprintf(stderr, "no distance function\n");
	free(pOffset);
	return NULL;
    }
    if (nbpoints <= LeafCapacity + 2){ /* leaf */
	off_t offset_mask = m->leaf_pgsize - 1;
	off_t page_mask = ~(m->leaf_pgsize - 1);

	uint8_t ntype = 0;
	MVPFile m2;
	ph_mvp_init(&m2);

	/* open new file */
	char extfile[256];
	sprintf(extfile, "%s%d.mvp", m->filename, m->nbdbfiles);
	
	m2.fd = open(extfile, O_CREAT|O_RDWR|O_APPEND, 00777);
	if (m2.fd < 0){
	    free(pOffset);
	    return NULL;
	}
        m2.hash_type = m->hash_type;

	struct stat fileinfo;
	fstat(m2.fd, &fileinfo);

	/* test file size and open new file if necessary */
	if (fileinfo.st_size >= MaxFileSize){
	    if (close(m2.fd) < 0){
		perror("fclose");
	    }
	    m->nbdbfiles++;
	    sprintf(extfile, "%s%d.mvp", m->filename, m->nbdbfiles);
	    
	    m2.fd = open(extfile,O_CREAT|O_RDWR|O_APPEND, 00777);
	    if (m2.fd < 0){
		perror("open");
		free(pOffset);
		return NULL;
	    }
	}

	m2.file_pos = lseek(m2.fd, 0, SEEK_END);

	if (ftruncate(m2.fd, m2.file_pos + m->leaf_pgsize) < 0){
	    perror("ftruncate");
	}
	
	off_t end_pos = m2.file_pos + m->leaf_pgsize;
	pOffset->fileno = m->nbdbfiles;
	pOffset->offset = m2.file_pos;

	off_t pa_offset = m2.file_pos & page_mask;
        m2.buf = (char*)mmap(NULL,m->leaf_pgsize,
                                 PROT_READ|PROT_WRITE,MAP_SHARED,m2.fd,pa_offset);
	if (m2.buf == MAP_FAILED){
	    perror("mmap");
	    close(m2.fd);
	    free(pOffset);
	    return NULL;
	}

	m2.buf[m2.file_pos & offset_mask] = ntype;
	m2.file_pos++;

	/* find vantage points, sv1 and sv2 */
	DP *sv1 = points[0];
	DP *sv2 = NULL;
	float d, max_dist = 0;
	int max_pos = 0;
	for (int i=1; i< nbpoints;i++){
	    d = hashdist(sv1, points[i]);
	    if (d > max_dist){
		max_dist = d;
		max_pos = i;
	    }
	}
	
	if (max_pos > 0){
	    sv2 = points[max_pos]; /* sv2 is furthest point from sv1 */
	}

	ph_save_datapoint(sv1, &m2);
	ph_save_datapoint(sv2, &m2);

	m2.buf[m2.file_pos & offset_mask] = (uint8_t)Np;
	m2.file_pos++;

	/* write the leaf points */
	float d1, d2;
	off_t curr_pos = m2.file_pos;
	off_t last_pos = curr_pos + LeafCapacity*(2*sizeof(float) + sizeof(off_t));
	off_t dp_pos;
	for (int i=1;i<nbpoints;i++){
	    if (i==max_pos) /* skip sv2 */
		continue;
	    /* write d1, d2 */
	    d1 = hashdist(sv1, points[i]);
	    d2 = hashdist(sv2, points[i]);
            memcpy(&(m2.buf[m2.file_pos & offset_mask]), &d1, sizeof(float));
	    m2.file_pos += sizeof(float);
            memcpy(&(m2.buf[m2.file_pos & offset_mask]), &d2, sizeof(float));
	    m2.file_pos += sizeof(float);
	    
	    /* write the point[i] at the end and return to write what the offset is*/
	    curr_pos = m2.file_pos;
	    m2.file_pos = last_pos;
	    
	    dp_pos = ph_save_datapoint(points[i], &m2);
	    last_pos = m2.file_pos;
	    m2.file_pos = curr_pos;

	    memcpy(&(m2.buf[m2.file_pos & offset_mask]), &dp_pos, sizeof(off_t));
	    m2.file_pos += sizeof(off_t);
	}

	if (msync(m2.buf, sysconf(_SC_PAGE_SIZE), MS_SYNC) < 0){
	    perror("msync");
	}
	if (munmap(m2.buf, sysconf(_SC_PAGE_SIZE)) < 0)
	    perror("munmap");

	if (ftruncate(m2.fd, end_pos) < 0)
	    perror("ftruncate");

	if (close(m2.fd) < 0)
	    perror("fclose");
    } else {
	off_t offset_mask = m->internal_pgsize - 1;
	off_t page_mask = ~(m->internal_pgsize - 1);

	pOffset->fileno = 0;
	pOffset->offset = m->file_pos;

	off_t orig_pos = m->file_pos;
        off_t end_pos;

        if ((level > 0) && (saveall_flag == 1)){ /* append new page to mainfile, unmap/map to it */

	    if (msync(m->buf, m->internal_pgsize ,MS_SYNC) < 0){
		perror("msync");
		return NULL;
	    }
	    
	    if (munmap(m->buf, m->internal_pgsize) < 0){
		perror("munmap");
	    }

	    m->file_pos = lseek(m->fd, 0, SEEK_END);
            end_pos = m->file_pos + m->internal_pgsize;
	    if (ftruncate(m->fd, m->file_pos + m->internal_pgsize) < 0){
		perror("ftruncate");
		return NULL;
	    }
	    off_t pa_offset = m->file_pos & page_mask;
	    m->buf = (char*)mmap(NULL,m->internal_pgsize,PROT_READ|PROT_WRITE,MAP_SHARED,m->fd,pa_offset);
	    if (m->buf == MAP_FAILED){
		perror("mmap");
		free(pOffset);
		return NULL;
	    }
	    pOffset->offset = m->file_pos;
	}

	uint8_t ntype = 1;
	memcpy(&m->buf[m->file_pos++ & offset_mask],  &ntype, sizeof(uint8_t));
	
	/* choose vantage points, sv1, sv2 */
	DP *sv1 = points[0]; 
	DP *sv2 = NULL;
	float max_dist = 0.0, min_dist = 1000.0;
	float *dist = (float*)malloc(nbpoints*sizeof(float));
	if (!dist){
	    free(pOffset);
	    return NULL;
	}
	int max_pos = 0;
	for (int i=0;i<nbpoints;i++){
	    dist[i] = hashdist(sv1, points[i]);
	    if (dist[i] > max_dist){
		max_pos = i;
		max_dist = dist[i];
	    }
	    if ((dist[i] < min_dist) && (dist[i] != 0)){
		min_dist = dist[i];
	    }
            if (level <= PathLength)
		points[i]->path[level] = dist[i];
	}
	sv2 = points[max_pos]; /* sv2 is furthest away from sv1 */

	/* save sv1, sv2 */
	ph_save_datapoint(sv1, m);
	ph_save_datapoint(sv2, m);

	/* 1st tier pivots, M1, derived from the distance of each point from sv1*/
	float step = (max_dist - min_dist)/BranchFactor;
        float incr = step;

	float *M1 = (float*)malloc(LengthM1*sizeof(float));
	float *M2 = (float*)malloc(LengthM2*sizeof(float));
	if (!M1 || !M2){
	    printf("unable to allocate M1[] or M2[]\n");
	    free(pOffset);
	    free(dist);
	    return NULL;
	}


	for (int i=0;i<LengthM1;i++){
	    M1[i] = min_dist + incr;
	    incr += step;
	    memcpy(&(m->buf[m->file_pos & offset_mask]),&M1[i], sizeof(float));
	    m->file_pos += sizeof(float);
	}

	/*set up 1st tier sorting bins - contain pointers to DP points[] param so that we only 
          move pointers, not the actual datapoints */
	DP ***bins = (DP***)malloc(BranchFactor*sizeof(DP***));
	if (!bins){
	    fprintf(stderr, "mem alloc error\n");
	    free(pOffset);
	    free(dist);
	    free(M1);
	    free(M2);
	    return NULL;
	}
	int *mlens = (int*)calloc(BranchFactor, sizeof(int)); /*no. points in each bin */
	if (!mlens){
	    free(pOffset);
	    free(dist);
	    free(M1);
	    free(M2);
	    free(bins);
	    return NULL;
	}

	for (int i=0;i<BranchFactor;i++){
	    bins[i] = (DP**)malloc(Np*sizeof(DP**)); /*Np should be more than enough */            
	    if (!bins[i]){
		free(pOffset);
		free(dist);
		free(M1);
		free(M2);
		free(bins);
		free(mlens);
		return NULL;
	    }
	}

	/* sort points into bins (except sv1 and sv2 )*/
	for (int i=1;i<nbpoints;i++){
	    if (i == max_pos)
		continue;
	    float cur_dist = dist[i];
	    /* check if <= M1[i] */
	    for (int j=0;j < LengthM1;j++){
		if (cur_dist <= M1[j]){
		    bins[j][mlens[j]] = points[i];
		    mlens[j]++;
		    break;
		}
	    }
	    /* check if > last M1[] pivot */
	    if (cur_dist > M1[LengthM1-1]){
		bins[BranchFactor-1][mlens[BranchFactor-1]] = points[i];
		mlens[BranchFactor-1]++;
	    }
	}

	/* print 1st level sort bins 
	for (int i=0; i<BranchFactor;i++){
	    int row_len = mlens[i];
	    for (int j=0;j<row_len;j++){
		printf(" %d %d %s\n", i, j, bins[i][j]->id);
	    }
	}
	*/
	
	/* set up 2nd tier sorting bins */
	/* each row from bins to be sorted into bins2, each in turn */
	DP ***bins2 = (DP***)malloc(BranchFactor*sizeof(DP***));
	if (!bins2){
	    free(pOffset);
	    free(dist);
	    free(M1);
	    free(M2);
	    free(bins);
	    free(mlens);
	    return NULL;
	}
	int *mlens2 = (int*)calloc(BranchFactor, sizeof(int)); /*number points in each bin */
	if (!mlens2){
	    free(pOffset);
	    free(dist);
	    free(M1);
	    free(M2);
	    free(bins);
	    free(bins2);
	    free(mlens);
	    return NULL;
	}

	for (int i=0;i<BranchFactor;i++){
	    bins2[i] = (DP**)malloc(Np*sizeof(DP**)); /* Np is more than enough */
	    if (!bins2[i]){
		free(pOffset);
		free(dist);
		free(M1);
		free(M2);
		free(bins);
		free(mlens);
		free(bins2);
		free(mlens2);
		return NULL;
	    }
	}
	
	
	off_t m2_pos = m->file_pos; /* start of M2 pivots */
	off_t child_pos = m2_pos + LengthM2*sizeof(float); /*pos where child offsets are written*/
	off_t last_pos = child_pos + Fanout*(sizeof(uint8_t) + sizeof(off_t)); /* last pos in
                                                                      in internal node */

	/* for each row of bin, sort the row into bins2 */
	for (int i=0;i < BranchFactor;i++){
	    int row_len = mlens[i]; /* length of current row, bins[i] */
	    for (int j=0;j < BranchFactor;j++){ /* reset the lengths to 0 */
		mlens2[j] = 0;
	    }
	    float *dist2 = (float*)realloc(dist,row_len*sizeof(float));
	    if (!dist2){
		free(pOffset);
		free(dist);
		free(M1);
		free(M2);
		free(bins);
		free(mlens);
		free(bins2);
		free(mlens2);
		return NULL;
	    }
	    dist = dist2;

	    /* 2nd tier pivots M2[], for row */
	    max_dist = 0;
	    min_dist = 1000;
	    for (int j=0;j<row_len;j++){
		dist[j] = hashdist(sv2, bins[i][j]);
		if ( dist[j] > max_dist)
		    max_dist = dist[j];
		if (dist[j] < min_dist)
		    min_dist = dist[j];
		if (level < PathLength){
		    bins[i][j]->path[level+1] = dist[j];
		}
	    }

	    step = (max_dist - min_dist)/BranchFactor;
	    incr = step;
	    
	    for (int j=0;j < LengthM1;j++){
		M2[j+i*LengthM1] = min_dist + incr;
		memcpy(&(m->buf[m2_pos & offset_mask]),&M2[j+i*LengthM1],sizeof(float));
		incr += step;
		m2_pos += sizeof(float);
	    }

	    /* sort bins[i] into bins2 */
	    for (int j=0;j < row_len; j++){
		DP *current = bins[i][j];
		/*check <= each M2 pivot  */
		for (int k=0;k<LengthM1;k++){
		    if (dist[j] <= M2[k+i*LengthM1]){
			bins2[k][mlens2[k]] = current;
			mlens2[k]++;
			break;
		    }
		}
		/* check > last M2 pivot  */
		if (dist[j] > M2[LengthM1-1 + i*LengthM1]){
		    bins2[BranchFactor-1][mlens2[BranchFactor-1]] = current;
		    mlens2[BranchFactor-1]++;
		}
	    }
	    
	    /* print 2nd tier sort bins 
	    for (int j=0;j<BranchFactor;j++){
		int r2 = mlens2[j];
		for (int k=0;k<r2;k++){
		    printf(" %d %d %d %s\n", i, j, k, bins2[j][k]->id);
		}
	    }
	    */
	    /* save child nodes */
	    FileIndex *pChild = NULL;
	    for (int j=0;j<BranchFactor;j++){
		m->file_pos = last_pos;
		pChild = ph_save_mvptree(m, bins2[j], mlens2[j], saveall_flag, level+2);
		if (pChild){ /* write filenumber and offset of child node */
		    last_pos = m->file_pos;
		    m->file_pos = child_pos;
		    memcpy(&(m->buf[m->file_pos & offset_mask]),&(pChild->fileno),sizeof(uint8_t));
		    m->file_pos++;
		    memcpy(&(m->buf[m->file_pos & offset_mask]), &(pChild->offset),sizeof(off_t));
		    m->file_pos += sizeof(off_t);
		    child_pos = m->file_pos;
		    
		} else { /* child node is null */
		    last_pos = m->file_pos;
		    m->file_pos = child_pos;
		    uint8_t emptyfileno = 0;
		    memcpy(&(m->buf[m->file_pos & offset_mask]), &emptyfileno, sizeof(uint8_t));
		    m->file_pos++;
		    off_t emptyoffset = 0;
		    memcpy(&(m->buf[m->file_pos & offset_mask]), &emptyoffset, sizeof(off_t)); 
		    m->file_pos += sizeof(off_t);
		    child_pos = m->file_pos;
		}
	    }
	    m->file_pos = last_pos;
	}

        /* remap to orig_pos */
	if ((level > 0) && (saveall_flag == 1)){
            /* unmap/remap to page with original position */
	    if (msync(m->buf, m->internal_pgsize, MS_SYNC) < 0)
		perror("msync");

	    if (munmap(m->buf, m->internal_pgsize) < 0)
		perror("munmap");

	    m->buf=(char*)mmap(NULL,m->internal_pgsize,PROT_WRITE,MAP_SHARED,m->fd,orig_pos & page_mask);
	    if (m->buf == MAP_FAILED){
		perror("mmap");
	    }
	    m->file_pos = orig_pos;
	}
	/* cleanup */
	free(bins);
	free(bins2);
	free(mlens);
	free(mlens2);
	free(dist);
	free(M1);
	free(M2);
    }
    return pOffset;
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

MVPRetCode ph_add_mvptree(MVPFile *m, DP *new_dp, int level){

    uint8_t ntype;
    off_t offset_mask, page_mask;
    hash_compareCB hashdist = m->hashdist;
    if (m->isleaf){
	offset_mask = m->leaf_pgsize - 1;
	page_mask = ~(m->leaf_pgsize - 1);
    } else {
	offset_mask = m->internal_pgsize - 1;
	page_mask = ~(m->internal_pgsize - 1);
    }
    off_t start_pos = m->file_pos;

    printf("add: level = %d, filepos = %ld, fd = %d\n", level,start_pos,m->fd);
    memcpy(&ntype, &m->buf[m->file_pos++ & offset_mask], sizeof(uint8_t));
    printf("ntype = %u\n", ntype);
    if (ntype == 0){
	uint8_t Np = 0;
	printf("add: leaf\n");
	DP *sv1 = ph_read_datapoint(m);
	if (sv1){
	    DP *sv2 = ph_read_datapoint(m);
	    if (sv2){
		off_t Np_pos = m->file_pos;
		memcpy(&Np,&m->buf[m->file_pos & offset_mask], sizeof(uint8_t));
		m->file_pos++;
 
		float d1 = hashdist(sv1,new_dp);
		float d2 = hashdist(sv2,new_dp);
		
		printf("sv1: %s, d1 = %f\n", sv1->id, d1);
		printf("sv2: %s, d2 = %f\n", sv2->id, d2);

		printf("Np = %u, at Np_pos = %ld\n", Np, Np_pos);

		off_t curr_pos, new_pos, point_pos;
		if (Np == 0){
                    printf("adding first point to leaf\n");
		    
		    memcpy(&m->buf[m->file_pos & offset_mask], &d1, sizeof(float));
		    m->file_pos += sizeof(float);
		    memcpy(&m->buf[m->file_pos & offset_mask], &d2, sizeof(float));
		    m->file_pos += sizeof(float);

		    curr_pos = m->file_pos;

		    m->file_pos += sizeof(off_t);
		    m->file_pos += (m->leafcapacity-1)*(2*sizeof(float) + sizeof(off_t));
		    new_pos = ph_save_datapoint(new_dp, m);
		    
		    memcpy(&m->buf[m->file_pos & offset_mask], &new_pos, sizeof(off_t));
		    
		    Np = 1;
		    printf("Np = %u at Np_pos = %ld\n", Np,Np_pos);
		    memcpy(&m->buf[Np_pos & offset_mask], &Np, sizeof(uint8_t));
		    
		} else if (Np < m->leafcapacity){
		    printf("adding point to leaf\n");
		    m->file_pos += (Np-1)*(2*sizeof(float)+ sizeof(off_t));
		    m->file_pos += 2*sizeof(float);
		    memcpy(&point_pos,&m->buf[m->file_pos & offset_mask], sizeof(off_t));
		    m->file_pos += sizeof(off_t);
		    
		    curr_pos = m->file_pos;
		    m->file_pos = point_pos;
		    uint8_t active;
		    uint16_t byte_len;
		    memcpy(&active,&m->buf[m->file_pos & offset_mask], sizeof(uint8_t));
		    m->file_pos++;
		    memcpy(&byte_len,&m->buf[m->file_pos & offset_mask], sizeof(uint16_t));
		    m->file_pos += sizeof(uint16_t);
		    printf("active = %u, byte length = %u\n", active, byte_len);
		    m->file_pos += byte_len;

		    new_pos = ph_save_datapoint(new_dp, m);
		    memcpy(&m->buf[curr_pos & offset_mask],&d1, sizeof(float));
		    curr_pos += sizeof(float);
		    memcpy(&m->buf[curr_pos & offset_mask],&d2, sizeof(float));
		    curr_pos += sizeof(float);
		    memcpy(&m->buf[curr_pos & offset_mask],&new_pos, sizeof(off_t));
		    curr_pos += sizeof(off_t);
		    
		    Np++;
		    printf("Np = %u, Np_pos = %ld\n", Np, Np_pos);
		    memcpy(&m->buf[Np_pos & offset_mask], &Np, sizeof(uint8_t));
		} else {
		    printf("Np now exceeds leaf capacity, creating new internal node w/nbpoints = %d\n",Np+3);
		    DP **points = (DP**)malloc((Np+3)*sizeof(DP**));
		    points[0] = sv1;
		    points[1] = sv2;
		    points[2] = new_dp;
		    
		    for (int i=0;i<Np;i++){
			m->file_pos += 2*sizeof(float);
			memcpy(&point_pos,&m->buf[m->file_pos & offset_mask], sizeof(off_t));
			m->file_pos += sizeof(off_t);
			
			curr_pos = m->file_pos;
			m->file_pos = point_pos;
			points[i+3] = ph_read_datapoint(m);
			m->file_pos = curr_pos;

			printf("point[%d]: %s\n", i+3, points[i+3]->id);
		    }
		    m->file_pos = start_pos;
		    if (!ph_save_mvptree(m, points, Np+3, 0, level+2)){
			fprintf(stderr, "unable to save new node\n");
		    }

		}
	    } else { /* put new point into sv2 pos */
	        printf("sv2 is null, add as sv2 for leaf, file_pos = %ld\n",m->file_pos);
		m->file_pos = start_pos;
		ntype = 1;
		Np = 0;


	    }
	} 
    } else if (ntype == 1){
	printf("internal\n");

	int LengthM1 = m->branchfactor - 1;
	int LengthM2 = (m->branchfactor)*LengthM1;

	DP *sv1 = ph_read_datapoint(m);
	DP *sv2 = ph_read_datapoint(m);
	
	float d1 = hashdist(sv1, new_dp);
	float d2 = hashdist(sv2, new_dp);

	printf("sv1: %s %f\n", sv1->id, d1);
	printf("sv2: %s %f\n", sv2->id, d2);

	float *M1 = (float*)malloc(LengthM1*sizeof(float));
	if (!M1){
	    fprintf(stderr,"mem alloc of M1[]\n");
	    return PH_MEMALLOC;
	}
	float *M2 = (float*)malloc(LengthM2*sizeof(float));
	if (!M2){
	    fprintf(stderr,"mem alloc of M2[]\n");
	    return PH_MEMALLOC;
	}

	memcpy(M1, &m->buf[m->file_pos & offset_mask], LengthM1*sizeof(float));
	m->file_pos += LengthM1*sizeof(float);

	memcpy(M2, &m->buf[m->file_pos & offset_mask], LengthM2*sizeof(float));
	m->file_pos += LengthM2*sizeof(float);

	for (int i=0;i<LengthM1;i++){
	    printf("M1[%d]=%f\n", i, M1[i]);
	}
	for (int i=0;i<LengthM2;i++){
	    printf("M2[%d]=%f\n", i,M2[i]);
	}


	if (level < m->pathlength)
	    new_dp->path[level] = d1;
	if (level < m->pathlength - 1)
	    new_dp->path[level+1] = d2;

	int pivot1, pivot2;
	uint8_t filenumber;
	off_t child_pos, curr_pos;
	off_t start_pos = m->file_pos;
	off_t orig_pos;
	
	/* check <= each M1 pivot */
	for (pivot1=0;pivot1 < LengthM1;pivot1++){
	    if (d1 <= M1[pivot1]){
		printf("d1 <= M1[%d]\n", pivot1);

		/* check <= each M2 pivot */
		for (pivot2 = 0; pivot2 < LengthM1;pivot2++){
		    if (d2 <= M2[pivot2+pivot1*LengthM1]){
			printf("d2 <= M2[%d]\n",pivot2+pivot1*LengthM1);
			printf("child[%d]\n", pivot2+pivot1*m->branchfactor);

			/* determine pos from which to read filenumber and offset */
			curr_pos = start_pos + (pivot2+pivot1*m->branchfactor)*(sizeof(uint8_t)+sizeof(off_t));
			m->file_pos = curr_pos;
			memcpy(&filenumber,&m->buf[m->file_pos & offset_mask],sizeof(uint8_t));
			m->file_pos++;
			memcpy(&child_pos,&m->buf[m->file_pos & offset_mask], sizeof(off_t));
			m->file_pos += sizeof(off_t);
			
			printf("filenum = %d, pos = %ld, orig_pos = %ld,fd=%d\n", filenumber, child_pos,m->file_pos,m->fd);
			/* save position and remap to new file/position */
			orig_pos = m->file_pos;
			printf("orig_pos = %ld, fd=%d\n", orig_pos, m->fd);
			MVPFile *m2 = _ph_map_mvpfile(filenumber,child_pos,m);
			if (m2){
			    if (ph_add_mvptree(m2, new_dp, level+2) != 0){
				return PH_NOSAVEMVP;
			    }
			}
			printf("back to orig_pos = %ld, fd=%d\n", orig_pos,m->fd);
			_ph_unmap_mvpfile(filenumber, orig_pos, m, m2);
			
		    }
		}
		/* check > last M2 pivot */
		if (d2 > M2[LengthM1-1+pivot1*LengthM1]){
		    printf("d2 > M2[%d]\n", LengthM1-1+pivot1*LengthM1);
		    printf("child[%d]\n", m->branchfactor-1+pivot1*m->branchfactor);
		    curr_pos = start_pos + (m->branchfactor-1+pivot1*m->branchfactor)*(sizeof(uint8_t)+sizeof(off_t));
		    m->file_pos = curr_pos;
		    memcpy(&filenumber,&m->buf[m->file_pos & offset_mask], sizeof(uint8_t));
		    m->file_pos++;
		    memcpy(&child_pos,&m->buf[m->file_pos & offset_mask], sizeof(off_t));
		    m->file_pos += sizeof(off_t);
		    
		    printf("filenum = %d, pos = %ld, orig_pos = %ld,fd=%d\n", filenumber,child_pos,m->file_pos,m->fd);
		    orig_pos = m->file_pos;
		    printf("orig_pos = %ld\n", orig_pos);
		    MVPFile *m2 = _ph_map_mvpfile(filenumber, child_pos, m);
		    if (m2){
			if (ph_add_mvptree(m2, new_dp, level+2) != 0){
			    return PH_NOSAVEMVP;
			}
		    }
		    printf("back to orig_pos = %ld, fd=%d\n", orig_pos,m->fd);
		    _ph_unmap_mvpfile(filenumber,orig_pos, m, m2);
		}

	    }
	}
	/* check > last M1 pivot */
	if (d1 > M1[LengthM1-1]){
	    printf("d1 > M1[%d]\n", LengthM1);
	    /*check <= each M2 pivot */
	    for (pivot2=0;pivot2 < LengthM1; pivot2++){
		if (d2 <= M2[pivot2+LengthM1*LengthM1]){
		    printf("d2 <= M2[%d]\n", pivot2+LengthM1*LengthM1);
		    printf("child[%d]\n", pivot2+LengthM1*m->branchfactor);

		    curr_pos = start_pos + (pivot2+LengthM1*m->branchfactor)*(sizeof(uint8_t)+sizeof(off_t));
		    m->file_pos = curr_pos;
		    memcpy(&filenumber, &m->buf[m->file_pos & offset_mask], sizeof(uint8_t));
		    m->file_pos++;
		    memcpy(&child_pos,&m->buf[m->file_pos & offset_mask], sizeof(off_t));
		    m->file_pos += sizeof(off_t);
		    
		    printf("filenumber = %d, child_pos = %ld, orig pos = %ld,fd=%d\n", filenumber,child_pos, m->file_pos,m->fd);
		    orig_pos = m->file_pos;
		    printf("orig_pos = %ld,fd=%d\n", orig_pos,m->fd);
		    MVPFile *m2 = _ph_map_mvpfile(filenumber,child_pos,m);
		    if (m2){
			if (ph_add_mvptree(m2,new_dp,level+2) != 0){
			    return PH_NOSAVEMVP;
			}
		    }
		    printf("back to orig pos = %ld, fd=%d\n", orig_pos,m->fd);
		    _ph_unmap_mvpfile(filenumber, orig_pos, m, m2);
		}
	    }
	    
	    /* check > last M2 pivot */
	    if (d2 > M2[LengthM1-1+LengthM1*LengthM1]){
		printf("d2 > M2[%d]\n", LengthM1-1+LengthM1*LengthM1);
		printf("child[%d]\n", m->branchfactor-1+LengthM1*m->branchfactor);
		curr_pos = start_pos + (m->branchfactor- 1+LengthM1*m->branchfactor)*(sizeof(uint8_t) + sizeof(off_t));
		m->file_pos = curr_pos;
		memcpy(&filenumber, &m->buf[m->file_pos & offset_mask], sizeof(uint8_t));
		m->file_pos++;
		memcpy(&child_pos, &m->buf[m->file_pos & offset_mask], sizeof(off_t));
		m->file_pos += sizeof(off_t);
		
		printf("filenumber = %d, child_pos = %ld, orig_pos = %ld,fd=%d\n", filenumber, child_pos, m->file_pos,m->fd);
		orig_pos = m->file_pos;
		printf("orig_pos = %ld, fd=%d\n", orig_pos,m->fd);
		MVPFile *m2 = _ph_map_mvpfile(filenumber, child_pos, m);
		if (m2){
		    if (ph_add_mvptree(m2, new_dp, level+2) != 0)
			return PH_NOSAVEMVP;
		}
		printf("back to orig pos = %ld, fd=%d\n", orig_pos,m->fd);
		_ph_unmap_mvpfile(filenumber, orig_pos, m, m2);
	    }
	}
    } else {
	fprintf(stderr,"unknown node type %u\n",ntype);
    }
    return PH_SUCCESS;
}

int ph_add_mvptree(MVPFile *m, DP **points, int nbpoints){

    printf("add to mvpfile %s %d points\n", m->filename, nbpoints);

    if (m->internal_pgsize == 0)
	m->internal_pgsize = sysconf(_SC_PAGE_SIZE);

    if (m->leaf_pgsize == 0)
	m->leaf_pgsize = sysconf(_SC_PAGE_SIZE);

    /* check to see that the pg sizes are at least the size of host page size */
    off_t host_pgsize = sysconf(_SC_PAGE_SIZE);
    if ((m->internal_pgsize < host_pgsize) || (m->leaf_pgsize < host_pgsize)){
	return -1;
    }

    /* pg sizes must be a power of zero */
    if ((m->internal_pgsize) & (m->internal_pgsize - 1))
	return -1;
    if ((m->leaf_pgsize) & (m->leaf_pgsize - 1))
	return -1;
    

    /* open main file */
    char mainfile[256];
    sprintf(mainfile, "%s.mvp", m->filename);

    printf("opening %s \n", mainfile);

    m->fd = open(mainfile, O_RDWR);
    if (m->fd < 0){
	perror("open");
	return -1;
    }
    /* map to first page */
    m->file_pos = 0;
    m->buf  = (char*)mmap(NULL,m->internal_pgsize,PROT_READ|PROT_WRITE,MAP_SHARED,m->fd,m->file_pos);
    if (m->buf == MAP_FAILED){
	perror("mmap");
	return -1;
    }
    /* read header within first HeaderSize bytes */
    printf("read header\n");
    char tag[17];
    int version;;
    int int_pgsize;
    int leaf_pgsize;
    
    memcpy(tag, &m->buf[m->file_pos], 16);
    tag[16] = '\0';
    m->file_pos += 16;
    
    memcpy(&version, &m->buf[m->file_pos], sizeof(int));
    m->file_pos += sizeof(int);

    memcpy(&int_pgsize, &m->buf[m->file_pos], sizeof(int));
    m->file_pos += sizeof(int);

    memcpy(&leaf_pgsize, &m->buf[m->file_pos], sizeof(int));
    m->file_pos += sizeof(int);

    memcpy(&m->nbdbfiles, &m->buf[m->file_pos++], 1);

    memcpy(&m->branchfactor, &m->buf[m->file_pos++], 1);

    memcpy(&m->pathlength, &m->buf[m->file_pos++], 1);

    memcpy(&m->leafcapacity, &m->buf[m->file_pos++], 1);

    memcpy(&m->hash_type, &m->buf[m->file_pos++], 1);

    printf("tag: %s, version: %d, intpgsize = %d, leafpgsize = %d\n", 
                                     tag,version,int_pgsize,leaf_pgsize);
    printf("nbdbfiles = %u, bf = %u, pl = %u, leafcap = %u, hashtype = %u\n",
	   m->nbdbfiles,m->branchfactor, m->pathlength,m->leafcapacity, m->hash_type);

    m->isleaf = 0;/* first file is never a leaf */
    m->file_pos = HeaderSize;

    int nbsaved = 0;
    printf("adding %d points ...\n",nbpoints);
    for (int i=0;i<nbpoints;i++){
	printf(" %d saving: %s...\n",i,points[i]->id);
        m->file_pos = HeaderSize;
	if (ph_add_mvptree(m, points[i], 0) != 0){
	    fprintf(stderr, "unable to save point: %s\n", points[i]->id);
	}
	printf("nbsaved = %d\n", nbsaved);
	nbsaved++;
	printf("enter key to continue\n");
	getchar();
    }
    printf("saved %d out of %d points\n", nbsaved, nbpoints);

    if (msync(m->buf, m->internal_pgsize, MS_SYNC) < 0){
	perror("msync");
	return -1;
    }

    if (munmap(m->buf, m->internal_pgsize) < 0){
	perror("ph_save_mvptree");
    }
   
    if (close(m->fd) < 0){
	perror("close");
    }

    return nbsaved;

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
 
    const char *dir_name = argv[1];/* name of dir to retrieve image files */
    const char *filename = argv[2];/* name of file to save db */

    MVPFile mvpfile;
    ph_mvp_init(&mvpfile);
    mvpfile.filename = strdup(filename);
    mvpfile.hashdist = hammingdistance;
    mvpfile.hash_type = UINT64ARRAY;


    ulong64 tmphash = 0;
    
    int nbfiles = 0;
    printf("dir name: %s\n", dir_name);
    char **files = ph_readfilenames(dir_name,nbfiles);
    if (!files){
	printf("mem alloc error\n");
	exit(1);
    }
    printf("nbfiles = %d\n", nbfiles);
    DP **hashlist = (DP**)malloc(nbfiles*sizeof(DP*));
    if (!hashlist){
	printf("mem alloc error\n");
	exit(1);
    }

    int count = 0;
    for (int i=0;i<nbfiles;i++){
	printf("file[%d]: %s ", i, files[i]);
	if (ph_dct_imagehash(files[i],tmphash) < 0){
	    printf("could not get for file %s\n",files[i]);
	    continue;
	}
	printf("hash = %llx\n", tmphash);

        hashlist[count] = ph_malloc_datapoint(mvpfile.hash_type,mvpfile.pathlength);
        if (!hashlist[count]){
	    printf("mem alloc error\n");
	    exit(1);
	}
	hashlist[count]->id = strdup(files[i]);
        void *ptr_hash = malloc(8);
	if (!ptr_hash){
            printf("unable to allocate mem\n");
            exit(1);
        }
        hashlist[count]->hash = ptr_hash;
	ulong64 *ptr = (ulong64*)hashlist[count]->hash;
        *ptr = tmphash;
	hashlist[count]->hash_length = 1;
         count++;
    }

    printf("add files to file %s\n", filename);
    int n = ph_add_mvptree(&mvpfile, hashlist, count);
    printf("number saved %d out of %d\n", n,count);
    if (n <= 0){
	printf("unable to add points to %s\n", filename);
    }
for(int i =0; i < nbfiles; i++)
	ph_free_datapoint(hashlist[i]);

free(hashlist);

    return 0;
}
