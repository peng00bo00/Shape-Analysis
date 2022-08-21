# Problem 1
(a)

$$
\begin{aligned}
M_\mathbf{p} \mathbf{n}_\mathbf{p} &= \frac{1}{2 \pi} \int_{-\pi}^{\pi} \kappa_\theta \mathbf{t}_\theta \mathbf{t}_\theta^T \mathbf{n}_\mathbf{p} \ d \theta \\
&= \frac{1}{2 \pi} \int_{-\pi}^{\pi} \kappa_\theta \mathbf{t}_\theta \cdot 0 \ d \theta \\
&= \mathbf{0}
\end{aligned}
$$

So, the normal $\mathbf{n}_\mathbf{p}$ is an eigen vector of $M_\mathbf{p}$ with the corresponding eigen value of 0.

(b)

Suppose $\mathbf{t}_{\max}$ and $\mathbf{t}_{\min}$ are the direction of maximum and minimum curvature $\kappa_{\max}$ and $\kappa_{\min}$. Then the directional curvature could be expressed as:

$$
\kappa_\theta = \kappa_{\min} \cos^2 \theta + \kappa_{\max} \sin^2 \theta
$$

$$
\mathbf{t}_\theta = \mathbf{t}_{\min} \cos \theta + \mathbf{t}_{\max} \sin \theta
$$

Then, 

$$
\begin{aligned}
M_\mathbf{p} \mathbf{t}_{\min} &= \frac{1}{2 \pi} \int_{-\pi}^{\pi} \kappa_\theta \mathbf{t}_\theta \mathbf{t}_\theta^T \mathbf{t}_{\min} \ d \theta \\
&= \frac{1}{2 \pi} \int_{-\pi}^{\pi} (\kappa_{\min} \cos^2 \theta + \kappa_{\max} \sin^2 \theta) \cdot (\mathbf{t}_{\min} \cos \theta + \mathbf{t}_{\max} \sin \theta) \cdot \cos \theta \ d \theta \\
&= \frac{1}{2 \pi} \int_{-\pi}^{\pi} 
\kappa_{\min} \cos^4 \theta \ \mathbf{t}_{\min} +
\kappa_{\max} \cos^2 \theta \sin^2 \theta \ \mathbf{t}_{\min} +
\kappa_{\min} \cos^3 \theta \sin \theta \ \mathbf{t}_{\max} +
\kappa_{\max} \cos \theta \sin^3 \theta \ \mathbf{t}_{\max} \ d \theta \\
&= \frac{3}{8} \kappa_{\min} \mathbf{t}_{\min} + \frac{1}{8} \kappa_{\max} \mathbf{t}_{\min} + 0 + 0 \\
&= \bigg(\frac{3}{8} \kappa_{\min} + \frac{1}{8} \kappa_{\max} \bigg) \ \mathbf{t}_{\min}
\end{aligned}
$$

So, $\mathbf{t}_{\min}$ is an eigen vector of $M_\mathbf{p}$ with the corresponding eigen value of $\frac{3}{8} \kappa_{\min} + \frac{1}{8} \kappa_{\max}$. Similarly, $\mathbf{t}_{\max}$ is another eigen vector of $M_\mathbf{p}$ with the corresponding eigen value of $\frac{1}{8} \kappa_{\min} + \frac{3}{8} \kappa_{\max}$. 