packages <- c(
    "sf",
    "ggplot2",
    "rnaturalearth",
    "rnaturalearthdata",
    "ggspatial",
    "dplyr"
)

missing_packages <- packages[!packages %in% rownames(installed.packages())]
missing_packages
install.packages(missing_packages)

library(sf)
library(ggplot2)
library(rnaturalearth)
library(ggspatial)
library(dplyr)

populations <- data.frame(
    Population = c("NA", "NL", "VF", "C.Fin", "Upp", "E", "L", "Ka"),
    Locality = c(
        "Arvidsjaur", "Luleå", "Fredrika", "Vaasa",
        "Uppsala", "Ürjaste", "Pelči", "Kalmar"
    ),
    Country = c(
        "Sweden", "Sweden", "Sweden", "Finland",
        "Sweden", "Estonia", "Latvia", "Sweden"
    ),
    Region = c(
        "North Sweden", "North Sweden", "North Sweden", "East",
        "South Sweden", "East", "East", "South Sweden"
    ),
    Latitude_group = c(
        "North", "North", "North", "North",
        "South", "South", "South", "South"
    ),
    lon = c(
        19.166819, 22.154650, 18.397577, 21.615769,
        17.638889, 24.366408, 21.829725, 16.361631
    ),
    lat = c(
        65.590331, 65.584150, 64.073111, 63.096000,
        59.858819, 59.095161, 56.763692, 56.661569
    )
)

pop_sf <- st_as_sf(
    populations,
    coords = c("lon", "lat"),
    crs = 4326
)


world <- ne_countries(
    scale = "medium",
    returnclass = "sf"
)

map_area <- world %>%
    filter(admin %in% c(
        "Sweden", "Finland", "Norway", "Denmark",
        "Estonia", "Latvia", "Lithuania", "Russia"
    ))

p <- ggplot() +
    geom_sf(data = map_area, fill = "grey92", colour = "grey45", linewidth = 0.3) +
    geom_sf(
        data = pop_sf,
        aes(shape = Region, colour = Latitude_group),
        size = 3
    ) +
    geom_text(
        data = populations,
        aes(x = lon, y = lat, label = Population),
        nudge_y = 0.22,
        size = 3.2
    ) +
    coord_sf(
        xlim = c(10, 27),
        ylim = c(55.5, 66.8),
        expand = FALSE
    ) +
    annotation_scale(
        location = "bl",
        width_hint = 0.25
    ) +
    annotation_north_arrow(
        location = "tl",
        which_north = "true",
        style = north_arrow_fancy_orienteering
    ) +
    labs(
        x = "Longitude",
        y = "Latitude",
        shape = "Region"
    ) +
    scale_colour_manual(
        values = c(
            "North" = "#f35050", # light red
            "South" = "#3892ed" # light blue
        ),
        name = "Latitude"
    ) +
    theme_classic() +
    theme(
        legend.position = "right",
        axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
    )

ggsave(
    filename = "rana_arvalis_sampling_map.png",
    plot = p,
    width = 16,
    height = 14,
    units = "cm",
    dpi = 300
)

ggsave(
    filename = "rana_arvalis_sampling_map.pdf",
    plot = p,
    width = 16,
    height = 14,
    units = "cm"
)
