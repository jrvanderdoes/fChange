## R CMD check results

0 errors | 0 warnings | 1 note

* This is a re-submission. I have removed the paper in the description which
    is accepted but was not yet published so the doi was inactive. I corrected
    spelling of an author's name by removing accents.
* Previous submissions corrected the following. Testing runtimes reduced to 
    under 10 minutes. Added references to the description in addition to the 
    numerous ones in the package. Several general references are given to 
    prevent the description from becoming too long with pages of references. 
    Upon publication of other papers, this will be expanded. Added 
    documentation for functions that return nothing. Removed commented examples. 
    Return user preferences after modification using 'on.exit()'. 
* Additionally, I would like to thank the reviewers for some helpful links, 
    which directed me to make a few other minor changes in-line with CRAN policy.
* This is a new release.
* As discussed in a previous email, this is reconstruction to a previous posted 
    package under the same name. That package was archived on 2020-03-07 as 
    check problems were not corrected despite reminders. We likewise have been 
    unable to reach the original maintainer. However, one of the original 
    authors, Gregory Rice, is listed as an author on this package. The other, 
    Alexander Aue, has approved and is part of the GitHub repo even though he is 
    not an author of this package. Despite using the same name, the code is 
    entirely redesigned, treats data completely differently, and was written by
    the authors attached to this version. Reasons for using the same name are 
    outlined in the previous email, which can be forwarded as needed.
* Some tests revealed a note that Rcpp is not needed for 'LinkingTo' in the 
    description file. If this occurs, this is incorrect. Rcpp is necessary for 
    the current code. I will aim to convert the code by the next submission.
