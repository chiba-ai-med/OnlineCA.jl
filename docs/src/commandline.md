# OnlineCA.jl (Command line tool)

All functions can be performed as command line tool in shell window and same options in [OnlineCA.jl (Julia API)](@ref) are available.

After installation of `OnlineCA.jl`, command line tools are saved at `YOUR_HOME_DIR/.julia/v1.x/OnlineCA/bin/`.

The functions can be performed as below.

## Correspondence Analysis (CA)
```bash
shell> julia YOUR_HOME_DIR/.julia/v1.x/OnlineCA/bin/ca \
--input Data.zst \
--outdir OUTDIR \
--dim 3 \
--noversamples 5 \
--niter 3 \
--chunksize 100
```

## Sparse Correspondence Analysis (SparseCA)
```bash
shell> julia YOUR_HOME_DIR/.julia/v1.x/OnlineCA/bin/sparse_ca \
--input Data.mtx.zst \
--outdir OUTDIR \
--dim 3 \
--noversamples 5 \
--niter 3 \
--chunksize 100
```

## Binary COO Correspondence Analysis (BinCOOCA)
```bash
shell> julia YOUR_HOME_DIR/.julia/v1.x/OnlineCA/bin/bincoo_ca \
--input Data.bincoo.zst \
--outdir OUTDIR \
--dim 3 \
--noversamples 5 \
--niter 3 \
--chunksize 100
```
