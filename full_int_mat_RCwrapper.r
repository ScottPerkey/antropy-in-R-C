start_time<-Sys.time()
library(readr)

dyn.load("all_ent.so")

extract_column_tuple<- function(directory,column_name){
	file_list <- list.files( path= directory_path,pattern= "*.csv", full.names=TRUE  )
	first_15_files<-head(file_list,15)
	tuple_list<-list()
	for (file in first_15_files){
		data<-read.csv(file)
		if (column_name %in% colnames(data)){
			filtered_data<-data[!is.na(data[[column_name]])&data[[column_name]]!='',]
			#extract first entry from the first column
			if (nrow(filtered_data)>20) {
				first_entry<-data[1,1]
				#extract the entire column of the STRING column
				mag_column<- data[[column_name]]
				#create tuple of the form ((first_entry,coumn3))
				tuple<-list(first_entry,mag_column)
				tuple_list<-c(tuple_list,list(tuple))
			}
		}
	}
	return(tuple_list)
}

multent<-function(mag_list){
	#m = size of vectors so m = comparing 2
	#n is the length of the time series
	#r is the distance the standard is 0.2*stddev(list)
	#ME is the pass by reference of the final value
	y<-mag_list
	n<-length(mag_list )
	M<-2
	r<-0.2*sd(mag_list)
	ME<-as.double(c(0))
	max_scale<-6
	result_mult<-.C("MSE",as.double(y),as.integer(M),as.double(r),as.integer(n),as.integer(max_scale) , ME=ME)
#	print(result_mult$ME)
	return(result_mult$ME)
}

sampent<-function(mag_list ){
	#m = size of vectors so m = comparing 2
	#n is the length of the time series
	#r is the distance the standard is 0.2*stddev(list)
	#SampEn is the pass by reference of the final value
	y<-mag_list
	n<-length(mag_list )
	M<-2
	r<-0.2*sd(mag_list)
	SampEn<-as.double(c(0))
	result_samp<-.C("sampen",as.double(y),as.integer(M),as.double(r),as.integer(n),SampEn=SampEn )
	return(result_samp$SampEn)
}

fftw_transform<-function(input){
	.Call("fftw_transform", as.numeric(input))
}

fourier_abs_sorted<-function(list_1){
	fftw_reals<-fftw_transform(list_1)[,1]
	fftw_reals_abs<-abs(fftw_reals)
	sorted_fftw_abs<-sort(fftw_reals_abs, decreasing=TRUE)
	return(sorted_fftw_abs)
}

shuffle_list<-function(x){
	return_x<-.C("shuffle",x=as.double(x),as.integer(length(x)))
	return(return_x$x)
}

wasserstein_distance<-function(p,q){
	n<-as.integer(length(p))
	wass<-as.double(c(0))
	result_wass<-.C("wasserstein_distance_r",as.double(p), as.double(q), n, wass=wass)
	return(result_wass$wass )
}

fourier_wass_normalized<-function(mags_input_list){

	iterations<-2
	list_wass_case_1<-numeric(0)
	list_wass_case_2<-numeric(0)
	wass_list_four_norm<-numeric(0)
	for ( i in 1:iterations ){
		mags_permutate_1<-shuffle_list(mags_input_list)
		fft1<-fourier_abs_sorted(mags_input_list)
		#print(mags_permutate_1[[2]])
		fft2<-fourier_abs_sorted(mags_permutate_1)
		mags_permutate_2<-shuffle_list(mags_permutate_1)
		fft3<-fourier_abs_sorted(mags_permutate_2)
		list_wass_case_1<-c(list_wass_case_1,wasserstein_distance(fft1,fft2))
		list_wass_case_2<-c(list_wass_case_2,wasserstein_distance(fft2,fft3))
	}
	for (i in 1:length(list_wass_case_1)){
		wass_list_four_norm<-c(wass_list_four_norm,(list_wass_case_1/mean(list_wass_case_2)))
	}
	return(mean(wass_list_four_norm))
}
integration_entropy<-function( mags_input   ){
	x_values <-seq(0,length(mags_input)-1)
	fft1_4_int_y<-fourier_abs_sorted(mags_input)


	interp_function<-splinefun(x_values,fft1_4_int_y)

	#print(interp_function)
	full_int_ent<-integrate(interp_function, 0, length(fft1_4_int_y),subdivisions=1000)
	return(full_int_ent$value)
}

process_tuples<-function(full_mag_QSO){
	#-----------------------------------------------------------------
	#Magnitudes list of list

	mags<-full_mag_QSO[[2]]
	
	#-----------------------------------------------------------------
	#OID's 	

	OID <- full_mag_QSO[[1]]

	#-----------------------------------------------------------------
	#Sample Entropy in C
	sampent_c<-sampent( mags )
	#-----------------------------------------------------------------
	#Multi-scale Entropy in C

	multent_c<-multent(mags)
	print(multent_c)
	#-----------------------------------------------------------------
	#Fourier Coeff WD Normalized in C

	four_norm<-fourier_wass_normalized(mags)
	#print(four_norm)

	#-----------------------------------------------------------------
	#Integration fourier coefficient entropy 

	integ_int<-integration_entropy(mags)

	#-----------------------------------------------------------------
	#returning the full row

	return(list( OID,sampent_c,four_norm,integ_int ))
	#return(list( OID,sampent_c,multent_c,four_norm,integ_int ))
	#-----------------------------------------------------------------
}	






directory_path<-"/sharedata/fastdisk/anascc/LensedQuasars/alberto_QSO/RA_10.5635869145"
column_name<-"mag"
vec_list<-extract_column_tuple(directory_path,column_name)
fuck<-process_tuples(vec_list[[3]])
print(fuck)
#main<-function(){}
#main()
end_time<-Sys.time()
execution_time<-as.numeric(difftime( end_time, start_time, units="secs"  ))
cat("It took",execution_time,"seconds to run.\n")

