install.packages("tidyverse", dependencies = TRUE)
install.packages("sf", dependencies = TRUE)
install.packages("raster", dependencies = TRUE)
install.packages("spData", dependencies = TRUE)
if (!require("devtools")) install.packages("devtools")
devtools::install_github("Nowosad/spDataLarge")
install.packages("sp", dependencies = TRUE)


library(tidyverse)
library(sf)          # classes and functions for vector data
library(raster)      # classes and functions for raster data
library(spData)      # load geographic data
library(spDataLarge) # load larger geographic data
library(sp)
library(tmaptools)

#Map data prepaird
data(world)
rs7519753_world <- read.table("rs7519753_world.txt",sep = "\t",stringsAsFactors = F,header=T)


#Calculating country's center point
world_cents = st_centroid(world, of_largest = TRUE)

country_cents <- world_cents[(world_cents %>% st_set_geometry(NULL) %>% as.data.frame())[,"name_long"] 
                             %in% rs7519753_world$NAME,] %>% dplyr::select(name_long,geom)

#Fill in missing data
for (i in c(15:18)) {
  country_cents[i,1] <- rs7519753_world[i,1]
  country_cents$geom[i] <- st_point(c(rs7519753_world[i,4],rs7519753_world[i,5]))
}

country_cents <- country_cents[order(country_cents$name_long),]


#ggplot data preparation
origin_cols <- c("#0070C0","#00B050")

country_SNP <- rs7519753_world[rs7519753_world$NAME %in% as.data.frame(st_set_geometry(country_cents,NULL))[,"name_long"],]


country_SNP <- country_SNP[order(country_SNP$NAME),]

cents_data <- country_SNP %>% 
  mutate(FID= factor(1:n())) %>% 
  dplyr::select(NAME,FID,Proportion,ancestry,C,T) %>% 
  gather(key=origin, value=perc,C,T, factor_key=TRUE)



grobs_cents <- lapply(split(cents_data, cents_data$FID), function(x) {
  ggplotGrob(ggplot(x, aes(x="", y=-perc, fill=origin)) +
               geom_bar(width=1, stat="identity") +
               scale_y_continuous(expand=c(0,0)) +
               scale_fill_manual(values=origin_cols) +
               theme_ps(plot.axes = FALSE)+coord_polar(theta = 'y')+
               labs(title = paste(x$ancestry[1],x$Proportion[1],sep = ":"))+
               theme(plot.title = element_text(size = 8, colour = I("#E11524"))))
})


#Make ggplot plot and map plot have the same name
names(grobs_cents) <- country_cents$name_long


#Drawing
tm_shape(world)+
  tm_polygons(col = "#F0EEE3",border.col = NULL)+
  tm_shape(country_cents)+
  tm_symbols(size=1, shape="name_long",
             border.col = NULL,
             shapes=grobs_cents, 
             sizes.legend=c(.5, 1,3)*1e6, 
             scale=1, 
             legend.shape.show = FALSE, 
             legend.size.is.portrait = TRUE, 
             shapes.legend = 22, 
             title.size = "Population",
             group = "Charts",
             id = "name")