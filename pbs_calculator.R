# This script can take Fst Data from VCFTools
# and calculate PBS values

# Charging library
library("dplyr")
library("ggplot2")

## Read args from command line
args <- commandArgs(trailingOnly = T)

# Uncomment for debbuging
# Comment for production mode only
#args[1] <- "./dir_my_fst"
#args[2] <- "pbs_snp.png"
#args[3] <- "pbs.tsv"

## Place args into named object
file_dir <- args[1]
jpg_file <- args[2]
tsv_file <- args[3]

## Using function for reading
fst_reader_snp<-function(filename,pops){
  fst<-read.table(file = filename, header = T, stringsAsFactors = F, sep = "\t") #read in fst output from vcftools
  colnames(fst)[3]<-paste(pops,".Fst",sep = "") #change fst column name to identify pops
  fst[,3][which(is.na(fst[,3]))]<-0 #change NA's to 0, the NA's are produced by vcftools when there is no variation at a site
  assign(pops,fst,envir = .GlobalEnv) #export dataframe with pop names as variable names
}

## Reading data
snp.name <- list.files(path = file_dir, pattern = "*.fst", full.names = T)

## Preparing data
pop1_snp <- as.list(strsplit(snp.name[1], "/") [[1]])
pop1s_saved <- as.character(pop1_snp[[3]])

pop2_snp <- as.list(strsplit(snp.name[2], "/") [[1]])
pop2s_saved <- as.character(pop2_snp[[3]])

pop3_snp <- as.list(strsplit(snp.name[3], "/") [[1]])
pop3s_saved <- as.character(pop3_snp[[3]])

## Reading data
psnp1 <- fst_reader_snp(filename = snp.name[1], pops = pop1s_saved)
psnp2 <- fst_reader_snp(filename = snp.name[2], pops = pop2s_saved)
psnp3 <- fst_reader_snp(filename = snp.name[3], pops = pop3s_saved)

# Make a list of your dataframes and join them all together
predata <- left_join(x = psnp1,
                     y = psnp2,
                     by = c("CHROM", "POS"))

all_fst <- left_join(x = predata,
                     y = psnp3,
                     by = c("CHROM", "POS"))

## PBS function
dopbs<-function(pop1_out,pop1_pop2,pop2_out) {
  Tpop1_out= -log(1-pop1_out)
  Tpop1_pop2= -log(1-pop1_pop2)
  Tpop2_out= -log(1-pop2_out)
  pbs= (Tpop1_out + Tpop1_pop2 - Tpop2_out)/2
  pbs
}

# Running for my data
my_pbs <-dopbs(pop1_out = all_fst[4],
              pop1_pop2 = all_fst[3],
              pop2_out = all_fst[5])

# Turn PBS into data frame and merge
pbsresults<-all_fst[,1:2]
pbsresults<-cbind(pbsresults,my_pbs)
pbsresults <- pbsresults %>% 
  rename(PBS_value = paste0(pop2s_saved,".Fst"))

#set negative PBS values to 0 for convenience
pbsresults$PBS_value[which(pbsresults$PBS_value<0)]<-0

# Arranging data to see values
arreglado <- pbsresults[order(-pbsresults$PBS_value),]

## Saving df
write.table(arreglado, file = tsv_file, col.names = T, row.names = F, sep = "\t")

## Making manhattan plot
## vamos a quedarnos solo con datos que tengan cromosomas 1-22, X o Y
valid_chroms <- c(1:22, "X", "Y")

#limpiamos solo para quedarnos con cromosomas validos
clean_data.df <- arreglado %>%
  filter( CHROM %in% valid_chroms ) %>%
  ##Para poder ver chrom X y Y, los renumeramos a 23 y 24 respectivamente
  mutate( CHROM = gsub( pattern = "X", replacement = "23", x = CHROM ),
          CHROM = gsub( pattern = "Y", replacement = "24", x = CHROM )) %>%
  ## convertir en numericos los valores de chromosoma y posicion
  mutate( CHR_ID = as.numeric(CHROM),
          CHR_POS = as.numeric(POS))

##Preparamos la data para graficar en X un solo eje,
# a pesar de que sean diferentes cromosomas
single_axis.df <- clean_data.df %>% 
  # Compute chromosome size
  group_by( CHR_ID ) %>% 
  summarize(chr_len = max( CHR_POS) ) %>%
  # Calculate cumulative position of each chromosome
  mutate( chr_start_on_x = cumsum(chr_len) - chr_len) %>%
  select( -chr_len )

## Juntar clean data y el single axis para homogenizar la posicion
adjusted.df <- left_join( x = clean_data.df, y = single_axis.df,  by = "CHR_ID" ) %>%
  # ajustar posicion de acuerdo a un solo eje X
  arrange( CHR_ID , CHR_POS) %>%
  mutate( position_on_x = CHR_POS + chr_start_on_x,
          CHR_ID = as.factor(CHR_ID) )

## Vamos a crear un vector con las etiquetas para el eje X
xaxis.df = adjusted.df %>% group_by( CHR_ID ) %>%
  summarize( center = ( max( position_on_x ) + min( position_on_x) ) / 2 ) #Calculamos la posición central para cada cromosoma

##Limpiar dataframes que no se usaran
rm( clean_data.df )
rm( single_axis.df )

## Creating
man_uno.p <- ggplot( data = adjusted.df,
                     mapping = aes( x = position_on_x,
                                    y = PBS_value,     # transformamos los valores de Y a menos log10
                                    color = CHR_ID ) )   +
  geom_point( alpha = 0.2 )
man_uno.p

## Creating basic plot
# Vamos a cambiar la escala de colores, a algo mas estandar
man_dos.p <- man_uno.p +
  scale_color_manual( values = rep(c("#1B998B", "#E71D36"), 12 ))
man_dos.p

# vamos a poner los nombres de los cromosomas
man_tres.p <- man_dos.p +
  scale_x_continuous( label = xaxis.df$CHR_ID,      # Le ponemos un eje especial a X
                      breaks= xaxis.df$center ) 
man_tres.p

# Agregamos una linea roja con el corte de PBS
PBS_cutoff <- 0.2
man_cuatro.p <- man_tres.p +
  geom_hline( yintercept = PBS_cutoff,     
              color = "#2B2C28", lty = "dashed")        
man_cuatro.p

# Ponemos titulo, etiqueta en X y modificamos tema
man_cinco.p <- man_cuatro.p +
  labs(title = "PBS by SNP´s",    # ponemos titulo dinamico con la funcion paste
       x = "CHR",
       y = "PBS value") +
  theme_light(base_size = 12) +                              
  theme( legend.position="none",            
         panel.grid = element_blank(),
  )
man_cinco.p

# Separamos titulos de eje de plot
man_seis.p <- man_cinco.p +
  theme(axis.title.x = element_text(margin=margin(10,0,0,0), face = "bold", color = "grey20"),
        axis.title.y = element_text(margin=margin(0,10,0,0), face = "bold", color = "grey20"),
        plot.title=element_text(size=15,face="bold", color = "grey20"))
man_seis.p

# Guardar plot
ggsave( filename = jpg_file,     # nombre del archivo que se va a crear
        plot = man_seis.p,              # cual grafico guardamos?
        width = 15,                 # Ancho 10
        height = 8,                # Altura 10
        units = "in",               # Pulgadas, "in"ches
        dpi = 300 ) 
