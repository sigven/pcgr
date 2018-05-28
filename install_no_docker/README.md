### No-docker install

This is an alternative installation guide that does not require Docker:

```
bash -x install.sh
```

The script pulls dependencies from public repositories like conda and CRAN, and does not guarantee a successful install. Due to ongoing updates of the packages in such repositories, you might end up with version conflics, or unsupported versions for some packages, or missing packages for your system. Attempt on your own risk and try to stick to the dockerized version if possible.

After installing all dependencies, you will still need to download and uncompress data bundles, according to the [README](https://github.com/sigven/pcgr).