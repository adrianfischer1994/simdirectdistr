Instructions for the Autoencoder: 
1. Download the MNIST database and store it in a folder "data" in the same directory as "autoencoder.py".
2. Train the model with "autoencoder.py". It will save the model and the latent space representation for each digit in a directory "result".
3. Use "mnist_estimation_vMF.R" or "mnist_estimation_fisherbingham.R" to fit and generate from vMF or Fisher-Bingham. It will save the generated points on the sphere in a file "gen_sph.txt".
4. Use "num_save.py" in order to apply the decoder to the generated points on the sphere. It will save the results in a file "num.txt".
5. You can then use the Mathematica files to plot the generated handwritten digits or the latent representation on the sphere.

Parts of the code in "audoencoder.py" were taken from https://github.com/nicola-decao/s-vae-pytorch/, see [[1]](#citation)(http://arxiv.org/abs/1804.00891).

## References
```
[1] Davidson, T. R., Falorsi, L., De Cao, N., Kipf, T., and Tomczak, J. M. (2018).
Hyperspherical Variational Auto-Encoders.
34th Conference on Uncertainty in Artificial Intelligence (UAI-18).
```
