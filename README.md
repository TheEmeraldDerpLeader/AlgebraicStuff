# AlgebraicStuff
 Unfinished project for classifying real, unitial, associative algebras, and investigating weird ways to quotient their exponential groups

## Summary
Composed of two related projects. The first uses generalized Newton's method to find associative algebras and determine if two algebras are isomorphic via some linear transformation. Unfortunately hits a lot of false positives, e.g. literature says there are 6 unique 2D (not including unit) algebras, however this program will find more than 6 unless starting with the 6 algebras already assumed. I'm pretty sure this is a result of randomly generated algebras being very close to lower complexity algebras, which are more likely to fail the isomorphism test. E.g. the program randomly generates an alg which is equivalent to another under the transformation i => 1 + 0.0001i.

The second also uses generalized Newton's method to represent the exponential group of real unitial associative algebras, but under a novel quotient. The exponential group is normally defined via a summation exp(z) = sum z^n/n!. The new approach uses a similar summation, but instead recognizes a smaller set of vectors T called the translation space, with the orthogonal space R called the rotation space. For t an element of T and r an element of R, we have the group exp(t+r) = exp(t)\*exp(r), and we say that two elements of this group are the same if their ts are the same. Note that all elements of the group = some exp(t+r), even if the unit is removed from the rotation space. While convoluted and inefficient, the test cases on the quaternions and split-quaternions have shown that this method does in fact provide representations of S3 and H3, which contrast with the usual methods that involve the sandwich product and the typical exponential space. I have yet to do testing on other algebras to see how this holds up on non isotropic algebras

## Dependencies
The first project only depends on my [SymbolicCalculator](https://github.com/TheEmeraldDerpLeader/SymbolicCalculator) and [GLM](https://glm.g-truc.net/0.9.8/index.html). The second also depends on [SFML](https://www.sfml-dev.org). Can't be bothered to give more build instructions.

Not gonna be working more on this project anytime soon. It's a massive headache, both subprojects are just running into wall after wall. If I have some great idea on how to both separate false positives from good algebra candidates and reduce the computation time for these false positives, then I may come back to it.