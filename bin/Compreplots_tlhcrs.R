
#### script is designed to visualise the differences between TLHcrs repeats across a set of genomes

suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(ggtree)))
suppressMessages(suppressWarnings(library(phytools)))
suppressMessages(suppressWarnings(library(aplot)))
suppressMessages(suppressWarnings(library(dplyr)))


##create a simpler tree without the sokoriana root
raw_tree = ape::read.tree("assemblies.mashtree.bootstrap.dnd")

##root this one on but b y using the midroot function
midroot_tree=phytools::midpoint.root(raw_tree)

tree=suppressMessages(suppressWarnings(print(ggtree(midroot_tree)+
  geom_tiplab(as_ylab = TRUE, face="bold", align = T, size=0)+
  geom_nodelab(aes(label=label)))))

tip_order=tree$data %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label)


###plotting repeat similarities using global ANI in tile colours
within_gANI=read.csv(file="gANI.within_repeats.tsv", header=T, sep='\t')

within_gANI$gANI = within_gANI$average_gANI*100

##rearrange order of genomes by tree
#within_gANI$sample = factor(within_gANI$sample, levels = tip_order)

##rearrange order of genomes by tree
within_gANI$sample = factor(within_gANI$sample, levels = tip_order)
within=suppressMessages(suppressWarnings(print(ggplot(data=within_gANI, aes(x="", y=sample, fill= gANI)) +
  geom_tile(color = "black") +
  geom_text(aes(label = count), color = "white", fontface = "bold")+
  scale_fill_gradient2(high = "#FF0000", limits = c(0, 100))+
  scale_y_discrete(labels = setNames(within_gANI$sample, within_gANI$repeat_representative))+
  theme_classic()+
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), legend.position = "NA"))))

###
rep_gANI=read.csv(file="gANI.between_repeat_representatives.tsv", header=T, sep='\t')

rep_gANI$gANI = rep_gANI$gANI*100

##rearrange order of genomes by tree
rep_gANI$query_sample = factor(rep_gANI$query_sample, levels = tip_order)
rep_gANI$ref_sample = factor(rep_gANI$ref_sample, levels = tip_order)

between=suppressMessages(suppressWarnings(print(ggplot(data=rep_gANI, aes(x=query_sample, y=ref_sample, fill= gANI))+
  geom_tile(color = "black") +
  scale_fill_gradient2(high = "#FF0000", limits = c(0, 100))+
  theme_classic(base_size = 13) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "top", axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1),axis.text    = element_text(color = "black"), axis.ticks   = element_line(color = "black"), axis.line    = element_line(color = "black")))))


final=between %>% insert_left(tree, width=0.75) %>% insert_right(within, width=0.05)


##try to get dimensions of the plotted figure that adjust to the number of tips etc
ntips <- ape::Ntip(midroot_tree)
cell_size <- max(0.08, min(0.6, 15 / ntips))
tree_extra_width <- 4
##try to get a good width for the heatmap from within comparisons
side_extra_width <- ntips * cell_size 

# total figure dimensions
WGheight <- ntips * cell_size
WGwidth  <- (ntips * cell_size) + tree_extra_width + side_extra_width

ggsave("phylogeny_plus_gANI_heatmap.svg", plot = final, units = "in", height = WGheight, width = WGwidth, limitsize = FALSE)

