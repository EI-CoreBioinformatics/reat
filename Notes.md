# Release process

To create a new release, please use `bumpversion`. The following is an example of a release, when this command is executed it creates a commit and a tag which need to be pushed to the repository. 

```
bumpversion patch --verbose --tag --message "v0.5.1" --tag-message "v0.5.1

Docs update
Fixes issue with pre-processing files for EVM
Adds CLI parameter for pre-trained models for CodingQuarry, GlimmerHMM, and SNAP"
```

To push both the commit generated and the tag please use:

```
git push
git push --follow-tags
```

Once the tag and commit have been pushed a Github Action is triggered for the tag (based on the naming pattern vMajor.minor.patch) which will generate a release automatically.


