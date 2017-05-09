# unmarked_royle2004
Stan implementation of population density from unmarked surveys model (Royle 2004)

``validation.R`` -- script for simulating data and running models

``unmarked.stan`` -- basic unpooled model

``unmarked_pooled.stan`` -- hierachical version with beta pooling across detection rates ``p`` and exponential pooling across ``lambda`` latent rates

``unmarked_pooled2.stan`` -- hierarchical version with pooling across ``p`` and ``lambda``