# Problem 1
(a)

$$
\begin{aligned}
\frac{d}{dh} s[\gamma + h\mathbf{v}] \ \bigg|_{h=0} &= \frac{d}{dh} \int_0^1 \Vert \gamma' + h \mathbf{v}' \Vert_2 \ dt \ \bigg|_{h=0} \\
&= \int_0^1 \frac{d}{dh} \Vert  \gamma' + h \mathbf{v}' \Vert_2 \ dt \ \bigg|_{h=0} \\
&= \int_0^1 \frac{(\gamma' + h \mathbf{v}' ) \cdot \mathbf{v}'}{\Vert  \gamma' + h \mathbf{v}' \Vert_2} \ dt \ \bigg|_{h=0} \\
&= \int_0^1 \frac{\gamma'}{\Vert  \gamma'\Vert_2} \cdot \mathbf{v}' \ dt \\
&= \int_0^1 \mathbf{T}(t) \cdot \mathbf{v}'(t) \ dt
\end{aligned}
$$

(b)

$$
\begin{aligned}
\frac{d}{dh} s[\gamma + h\mathbf{v}] \ \bigg|_{h=0} &= \int_0^1 \mathbf{T}(t) \cdot \mathbf{v}'(t) \ dt \\
&= \langle \mathbf{T}, \mathbf{v}' \rangle \\
&\approx \nabla \gamma \cdot \mathbf{v}'
\end{aligned}
$$

(c)

When $\mathbf{v}(0) = \mathbf{v}(1) = \mathbf{0}$,

$$
\begin{aligned}
\frac{d}{dh} s[\gamma + h\mathbf{v}] \ \bigg|_{h=0} &= \int_0^1 \mathbf{T}(t) \cdot \mathbf{v}'(t) \ dt \\
&= [\mathbf{T} \cdot \mathbf{v}']_0^1 - \int_0^1 \mathbf{T}' \cdot \mathbf{v} \ dt \\
&= - \int_0^1 \mathbf{T}' \cdot \mathbf{v} \ dt \\
&= - \int_0^{s(1)} \mathbf{T}'(\bar{s}) \cdot \mathbf{v}(\bar{s}) \ d\bar{s} \\
&= - \int_0^{s(1)} \kappa(\bar{s}) \ \mathbf{N}(\bar{s}) \cdot \mathbf{v}(\bar{s}) \ d\bar{s}
\end{aligned}
$$

So, $\mathbf{w} = -\kappa(\bar{s}) \ \mathbf{N}(\bar{s})$.

# Problem 2
(a)

$$
s(\mathbf{x}) = \sum_{i=1}^{n-1} \Vert \mathbf{x}_{i+1} - \mathbf{x}_i \Vert_2
$$

(b)

$$
\begin{aligned}
\nabla_{\mathbf{x}_i} s &= \nabla_{\mathbf{x}_i} \bigg( \Vert \mathbf{x}_i - \mathbf{x}_{i-1} \Vert_2 + \Vert \mathbf{x}_{i+1} - \mathbf{x}_i \Vert_2 \bigg) \\
&= \nabla_{\mathbf{x}_i} \Vert \mathbf{x}_i - \mathbf{x}_{i-1} \Vert_2 + \nabla_{\mathbf{x}_i} \Vert \mathbf{x}_{i+1} - \mathbf{x}_i \Vert_2 \\
&= \frac{\mathbf{x}_i - \mathbf{x}_{i-1}}{\Vert \mathbf{x}_i - \mathbf{x}_{i-1} \Vert_2} + \frac{\mathbf{x}_i - \mathbf{x}_{i+1}}{\Vert \mathbf{x}_i - \mathbf{x}_{i+1} \Vert_2}
\end{aligned}
$$

So, the gradient of the length is the sum of (unit) directional vector of the two adjacent segments. Then, the norm of gradient is:

$$
\Vert \nabla_{\mathbf{x}_i} s \Vert_2 = 2 \sin \frac{\theta}{2}
$$

(d)

$$
\nabla_{\mathbf{x}_i} s = -2 \mathbf{N} \sin \frac{\theta}{2}
$$

$$
\begin{aligned}
\langle \nabla_{\mathbf{x}_i} s, \mathbf{v} \rangle &= -2 \sin \frac{\theta}{2} \ \mathbf{N} \cdot \mathbf{v} \\
&= \langle -\kappa \mathbf{N}, \mathbf{v} \rangle
\end{aligned}
$$

$$
\kappa = 2 \sin \frac{\theta}{2}
$$