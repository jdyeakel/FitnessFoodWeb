Methods
=======

We establish a framework that allows an individual forager to evaluate
trophic interactions over a great range of possibilities (species-rich
prey assemblages), where the dietary behavior is selected that maximizes
the individual’s current and future fitness. Our primary assumption that
underlies the mechanics of such decisions is that the forager evaluates
potential food items with respect to the similarities and differences
between a future potential resource, and the resource that it is
currently consuming. Accordingly, given the consumption of a specific
resource at time $t$, the preferred trophic interaction at time $t+1$
will vary over a range of prey-switching probabilities. If the preferred
resource at time $t+1$ is similar to the resource being consumed at time
$t$, the probability of switching is low. If the preferred resource at
time $t+1$ has qualities different than the resource being consumed at
time $t$, the probability of switching is high. By formulating the
dietary decision in this way, we focus on whether a forager maintains a
diet within a relatively small window of resource similarity, or whether
it frequently switches between different resources. Currently, we only
consider body size differences in potential resources to determine
‘similarity’, but this could potentially be expanded to include multiple
qualities.

Decision framework
------------------

The foraging decisions being weighted at time step $t$ ascribe unique
preference probabilities to the suite of available resources,
conditional on the forager’s current (focal) resource at time $t$,
$r_f(t)$. The decision framework is established on the probability that
a consumer preying upon a focal resource at time $t$ prefers to consume
a similar resource at time $t+1$ or switch to a different resource at
$t+1$. The suite of resources are thus defined in terms of their
similarity to each other, such that $\rm sim = 1$ denotes a resource
that is exactly similar to the focal resource, whereas $\rm sim = 0$
denotes a resource maximally different from the focal resource. As
similarity values range between 0 and 1, the preference probability for
the different resources can be thought of as the probability of
prey-switching, which varies according to a Beta Distribution. Thus,
preference probability distributions with lower expectations correspond
to higher probabilities of switching to a different resource, and
preference probability distributions with higher expectations correspond
to lower probability of switching to a different resource (i.e.
maintaining a diet similar to the current diet). The suite of foraging
decisions $D$ can thus be described by a range of beta distributions,
where the mean value increases from near 0 (high probability of
switching) to near 1 (low probability of switching); see Figure 1 for an
illustrative example.

State variables for single consumer scenario
--------------------------------------------

The state variables over which we determine the fitness-maximizing
foraging strategy include: 1) time $t$, 2) the forager’s current
energetic state $x(t)$, and 3) the current (focal) resource consumed by
the forager $r_f(t)$.

Derivation of the fitness function
----------------------------------

We will assume that fitness benefits increases with an individual’s
energetic state $x$, such that fitness at the terminal time interval
$t=T$ is some positive increasing function $\Phi(X(T))$. For time
periods $t < T$, fitness is maximized over the different decision
strategies $D$, where the choices of resources have probabilities
associated with different frequency distributions. The resource
preference frequency distribution assigns a preference probably to each
resource as a function of similarity to the focal resource
($\sigma_{r}={\rm sim}\{r_f,r\}$) is thus conditional on the decision
strategy $d \in D$ (illustrated in Figure 1) and the focal resource
$r_f$, thus denoted as $f(\sigma_{r}|d,r_f)$. The fitness increment due
to foraging $W$ is determined by the energetic state of the forager $x$,
as well as the energetic gains and losses associated with foraging for a
given resource. Fitness gain is thus a weighted average of potential
fitness benefits $W$ from consuming some resource and the probability of
finding $K=k$ items of the resource, which impacts the probability of
consumption. The probability that a consumer encounters a particular
quantity of a given resource has a Negative Binomial frequency
distribution $f(k|\xi_r,\nu)$, where $\xi_r$ is the mean encounter rate
of resource $r$, and $\nu$ is the dispersion parameter, which determines
‘patchiness’. We can thus incorporate the average encounter rate of a
resource, as well as the heterogeneity of the environment, in
determining the potential fitness consequences of foraging for a given
resource. It is initially assumed that the mean encounter rate declines
with increasing body size, and that dispersion is the same for each
resource (i.e. a product of the environment rather than a particular
species).

In addition to changes in fitness via foraging, we also allow each
species to obtain a reproductive fitness benefit in every time period
that it is alive (even if it dies within the given time period with
probability $m$), which is conditional on energetic state $\Psi(x)$. For
simplicity we will assume that the terminal fitness function is equal to
this fitness benefit over $x$ ($\Phi(x) = \Psi(x)$). Thus, our fitness
equation for time periods $t<T$ is written

$$\begin{aligned}
&F(x,r,t) = \\ \nonumber
&\underset{d}{\rm max}\sum_{p=r_{\rm min}}^{r_{\rm max}} f(\sigma_{p}|d,r_f) \left\{ \Psi(x) + (1 - m_r) \sum_{k=0}^{k_{\rm max}}f(k|\xi_p,\nu) W(x,p,t+1) \right\}.\end{aligned}$$

We assume that the probability of mortality in the resource encounter
$m_r$ increases with the body size of the resource, relative to
consumer’s body size, and that the maximum $m_r=0.5$. If $M_c$ and $M_r$
denote the body mass of the consumer and resource, and
$M_{\rm ratio}=M_r/M_c$, then

$$m_r = \frac{1}{2} - \left(\frac{1}{2}- m_r(M_{\rm ratio}=1)\right)^{M_{\rm ratio}^2},$$

where $m_r(M_{\rm ratio}=1)$ is a parameter for the probability of
mortality when the consumer and resource are of equal size. An organism
accrues energetic gains and costs as a function of its current energetic
state $x$, as well as the distribution and abundance of the resource in
the environment which influences how costly it is for the forager to
find and successfully consume a given resource. We assume that the
probability of a successful foraging bout $\rho$ increases with the
number of resources encountered, such that

$$\rho = 1 - \exp\left(-k^2/k_{\rm max}\right),$$

where $k$ is the number of resources encountered, and $k_{\rm max}$ is
the maximum possible number of resources encountered. We are now in a
position to define the energetic dynamics of foraging:

$$X(t+1) = X(t) + \rho(\mbox{gain-loss}) - (1-\rho)(\mbox{loss}),$$

where a successful foraging bout means that the consumer catches and
consumes a single resource (out of the $k$-sized group), and as the
group size increases, the probability of a successful capture likewise
increases. Furthermore, we might assume that the gains and losses are
allometrically derived and thus vary as a function of consumer and
resource body mass. In other words, gain may be of the form
$\epsilon M_r$ (where $\epsilon$ is transfer efficiency and $M_r$ is the
mass of the resource), while cost may be of the form $aM_c^b$, where $a$
and $b$ are allometric parameters and $M_c$ is the mass of the consumer.
The increment of fitness due to foraging is thus a function of energetic
state $x$, resource $r$, and time $t$, where

$$W(X(t)=x,r,t) = W \left( x + \rho(\epsilon M_r - aM_c^b) - (1-\rho)(aM_c^b),t+1 \right).$$

Initial Results
===============

Life history differences in diet-switching behaviors
----------------------------------------------------

(see Figure 2)

Early in life, before fitness has been accrued, the tendency to switch
increases with the body size of the focal prey. This behavior is
independent of the consumer’s energetic state, and primarily driven by
the increased mortality costs of larger prey, which eliminates any
future fitness gain. The progression from ’consistent’ to ’switching’ is
fast, occurring at intermediate resource body-size range.

Increasing the maximum probability of mortality does not alter the
decision matrix – it is still optimal to switch as the body sizes of
resources increase. In contrast, as the maximum probability of mortality
decreases, it is worthwhile to risk consumption of larger resources when
energetic state is low, as the potential gains are large relative to the
risk of mortality.

As fitness is accrued over an individual’s lifetime, the
fitness-maximizing behavior changes from avoiding resource-induced
mortality (to conserve future fitness gains) to increasing energetic
resources, which in turn leads to reproductive greater fitness gains per
time interval. Accordingly, when an individual nears the terminal time,
fitness-maximizing decisions become nonlinear relative to energetic
reserves and focal resource size. When the focal resource is
smaller-bodied, the likelihood of prey-switching decreases as energetic
reserves increase. When the focal resource is larger-bodied, the
likelihood of prey-switching increases as energetic resources increase.

Smaller animals have low mortality risk, but also low potential
energetic gain. Thus, if a consumer has a smaller-bodied focal resource
and its energetic reserves are high, reproductive fitness rewards are
already maximized, and switching to a larger resource only incurs
additional mortality risk. If it’s energetic reserves are low, it is
worth the additional mortality risk to switch to a resource that has a
greater chance of elevating the consumer’s energetic state and by proxy
it’s reproductive gain. Larger animals have a higher mortality risk and
large potential energetic gain. Accordingly, if the focal resource is
larger-bodied and the consumer’s energetic reserves are high,
reproductive fitness gains are already high and it is best to minimize
mortality risk by switching to a smaller resource. In contrast, if it’s
energetic resources are low, it is best to maintain it’s current focus
on larger-bodied resources with greater potential energetic gain (this
is true unless the focal resource is very large, such that mortality
risks are very great), which will ensure a higher probability of gaining
additional reproductive fitness gain before senescence.

These results are not strongly sensitive to the number of time-steps in
the model.

[These results are very tentative, as I’m not sure the fitness function
correctly deals with is issue]. When mortality is assumed to be
resource-independent, the dependence of foraging decisions on time
becomes minimal. In general, it is better to switch when consuming
smaller resources, and to maintain when consuming large resources. The
exception to this is when energetic resources become low: in this case,
the probability of switching moves from high to intermediate if the
focal resource body size is small, and energetic resources are low, due
to the higher probability of finding and catching smaller resources.

The effect of habitat heterogeneity
-----------------------------------

The qualitative fitness-maximizing decisions that emerge from the
fitness function are not strongly sensitive to habitat patchiness,
however there are quantitative differences in the decision solutions
between uniform vs. patchy environments. Uniform environments tend to
promote lower probabilities of diet-switching for smaller-sized focal
resources and for lower energetic states for the consumer. Heterogeneous
(patchy) environments tend to promote lower probabilities of
diet-switching for similarly low energetic states, but for larger-sized
focal resources.

It becomes better to switch from smaller-bodied prey to larger-bodied
prey when the environment is patchy because the potential rewards of
foraging for smaller prey are diminished as encounter rates become more
variable. In contrast, it becomes more beneficial to maintain diets on
larger-bodied resources when energetic reserves are low in patchy
environments, regardless of the greater mortality risks and lower
average encounter rates. Thus the window defining smaller probabilities
of diet-switching behaviors at lower energetic reserves shifts up the
focal resource body size axis. This indicates that specialization on
larger-bodied organisms provides greater fitness benefits in patchy
environments. This is due to the greater energetic gain attributable to
foraging on larger – albeit more rare and risky – body-sized resources.

Population-level impacts of diet-switching
------------------------------------------

From individual-level constraints to community structure
--------------------------------------------------------

[h!] ![image](Figure_Dist.pdf)

[fig~d~ist]

[h!] ![image](Figure_Patchy1.pdf)

[fig~p~atchy1]

