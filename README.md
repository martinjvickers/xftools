[![Build Status](https://travis-ci.org/martinjvickers/xftools.svg?branch=master)](https://travis-ci.org/martinjvickers/xftools)

# XFTOOLS

Tools written by Martin Vickers for use within the Xiaoqi Feng lab at the John Innes Centre.

## DISCLAIMER

Presently, these tools are still being tested and developed. There is little documentation outside of our lab for usage. In time as tools reach maturity I will add better documentation and better testing.

## Static Binary

For ease of use I have created static binaries which can be downloaded, extracted and run on Linux machines. These have been tested on Ubuntu 12.04 LTS, 14.04 LTS and 16.04 LTS as well as CentOS 6/7 and Scientific Linux 6/7. If these don't work on your machine you may need to compile from source. As I don't have a OSX machine I've not been able to test or create install instructions for OSX.

## Compile from source

TODO:

## Deployment

Deployments are managed automagically through Travis-CI. When we're ready for a new version ensure the following checklist is completed;

- [ ] New apps have a completed README.md, up-to-date version number and correction `--help` documentation/description
- [ ] Existing apps have the correct version number

After checklist, add the tag and push to origin

```
martin@x250:~/xftools$ git tag -a v0.0.4
martin@x250:~/xftools$ git push origin --tags
```

Post deployment checklist (Hopefully this can be automated)

- [ ] Deploy on the HPC
- [ ] Ensure that the static binaries work correctly on HPC
