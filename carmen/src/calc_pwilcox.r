write_to_max_p = function(max_val, max_p){
	for( i in (1:max_val) ){
		if( i < j ){
			for( j in ( i:(max_val)) ){
				legal = F
				for( test_val in (3: 10000) ){
					pval = pwilcox(test_val, i, j) * 2 # 2 tailed
					if( legal==F && pval < max_p ){
						legal = T
					}
					if( legal && pval > max_p ){
						write( c(i, j, test_val-1), 'c:\\code\\SingleQTL\\m_w_pvalues.txt', append=T )
						break
					}
				}
			}
		}
	}
	0
}


unlink('c:\\code\\SingleQTL\\m_w_pvalues.txt')

max_value = 40
max_p_value = 0.01
write( "-1 0.01 0", 'c:\\code\\SingleQTL\\m_w_pvalues.txt', append=T )

write_to_max_p(max_value, max_p_value)

write( "-1 0.001 0", 'c:\\code\\SingleQTL\\m_w_pvalues.txt', append=T )

max_p_value = 0.001
write_to_max_p(max_value, max_p_value)

write( "-1 0.0001 0", 'c:\\code\\SingleQTL\\m_w_pvalues.txt', append=T )

max_p_value = 0.0001
write_to_max_p(max_value, max_p_value)

write( "-1 0.00001 0", 'c:\\code\\SingleQTL\\m_w_pvalues.txt', append=T )

max_p_value = 0.00001
write_to_max_p(max_value, max_p_value)

write( "-1 0.000001 0", 'c:\\code\\SingleQTL\\m_w_pvalues.txt', append=T )

max_p_value = 0.000001
write_to_max_p(max_value, max_p_value)

write( "-1 0.0000001 0", 'c:\\code\\SingleQTL\\m_w_pvalues.txt', append=T )

max_p_value = 0.0000001
write_to_max_p(max_value, max_p_value)

write( "-1 0 0", 'c:\\code\\SingleQTL\\m_w_pvalues.txt', append=T )
