import sys

def extract_params(file):
	f = open(file, 'r')
	query = ""
	subject = ""
	index = 0
	
	for line in f:
		line = line.strip()

		if(index == 0):
			query = line 
		else:
			subject = line

		index += 1

	return query,subject

def split_query_subject(query, subject):
	l = len(query)
	m = len(subject)
	i = 0
	j = 0
	words_query = []
	words_subject = []
	
	#assumng the query and subject are of at least length 3
	while((i != l-1) and (i != l-2)):
		word_query = query[i] + query[i+1] + query[i+2]
		words_query.append(word_query)

		i += 1

	while((j != m-1) and (j != m-2)):
		word_subject = subject[j] + subject[j+1] + subject[j+2]
		words_subject.append(word_subject)

		j += 1

	return words_query,words_subject

def record_indices(query, subject):
	query_dict = {}
	subject_dict = {}
	index = 0

	for i in query:
		if(i in query_dict):
			query_dict[i].append(index) 
		else:
			query_dict[i] = [index]
		
		index += 1

	index = 0
	for i in subject:
		if(i in subject_dict):
			subject_dict[i].append(index)
		else:
			subject_dict[i] = [index]
		
		index += 1

	return query_dict,subject_dict

def align(query, subject, words_q, words_s, dict_q, dict_s, blosum_matrix, threshold):
	lenq = len(words_q)
	lens = len(words_s)
	hsp = {}
	drop = -3

	for windex in range(lenq):
		print(windex)
		word = words_q[windex]
		wr = windex + 3
		wl = windex - 1
		neighborhood_words = find_neighborhood_words(word, blosum_matrix, threshold)

		neighborhood_words = [word] + neighborhood_words

		for nword in neighborhood_words:
			
			if(nword in dict_s):
				indices = dict_s[nword]
				max_score = compute_score(word, nword, blosum_matrix)
				curr_score = max_score
				alignment = [word, nword]

				for i in indices:
					ir = i + 3
					il = i - 1

					while((wr < len(query)) and (ir < len(subject))):   #extend to the right

						if(str(query[wr] + subject[ir]) in blosum_matrix):
							pair = query[wr] + subject[ir]
						else:
							pair = subject[ir] + query[wr]

						pair_score = blosum_matrix[pair]

						if(pair_score < -3):
							break

						else:

							if((curr_score + pair_score) < (max_score + drop)):
								break

							else:
								curr_score += pair_score
								alignment[0] += query[wr]
								alignment[1] += subject[ir]

							if(curr_score >= max_score):
								max_score = curr_score
						wr += 1
						ir += 1
								

					while((wl >= 0) and (il >= 0)):
							
						if(str(query[wl] + subject[il]) in blosum_matrix):
							pair = query[wl] + subject[il]
						else:
							pair = subject[il] + query[wl]

						pair_score = blosum_matrix[pair]

						if(pair_score < -3):
							break

						else:
							curr_score += pair_score

							if(curr_score < (max_score + drop)):
								break

							else:
								alignment[0] = query[wl] + alignment[0]
								alignment[1] = subject[il] + alignment[1]

							if(curr_score >= max_score):
								max_score = curr_score

						wl -= 1
						il -= 1

					print('alignment= ' + str(alignment) + " " + str(curr_score))
					
					alignment = alignment[0] + "," + alignment[1]
					if(alignment not in hsp):
						hsp[alignment] = curr_score
					

	print(hsp)	
	return hsp			

def compute_score(word, word2, blosum_matrix):
	c0 = word[0] + word2[0]
	c1 = word[1] + word2[1]
	c2 = word[2] + word2[2]

	score0 = 0
	score1 = 0
	score2 = 0

	if(c0 in blosum_matrix):
		score0 = blosum_matrix[c0]

	else:
		score0 = blosum_matrix[word2[0] + word[0]]


	if(c1 in blosum_matrix):
		score1 = blosum_matrix[c1]

	else:
		score1 = blosum_matrix[word2[1] + word[1]]


	if(c2 in blosum_matrix):
		score2 = blosum_matrix[c2]

	else:
		score2 = blosum_matrix[word2[2] + word[2]]

	score = score0 + score1 + score2
	
	return score



def find_neighborhood_words(word, blosum_matrix, threshold):
	dict0 = [k for k,v in matrix_dict.items() if k.startswith(word[0]) or k.endswith(word[0])]
	dict1 = [k for k,v in matrix_dict.items() if k.startswith(word[1]) or k.endswith(word[1])]
	dict2 = [k for k,v in matrix_dict.items() if k.startswith(word[2]) or k.endswith(word[2])]
	word_list = []

	c0_match = word[0] + word[0]
	c0_score = blosum_matrix[c0_match]

	c1_match = word[1] + word[1]
	c1_score = blosum_matrix[c1_match]

	c2_match = word[2] + word[2]
	c2_score = blosum_matrix[c2_match]

	total_score = c1_score + c2_score
	for k in dict0:
		v = blosum_matrix[k]
		if(k != c0_match and ((total_score + v) >= threshold)):
			if(word[0] == k[0]):
				k = k[1]
			else:
				k = k[0]
			w = k + word[1] + word[2]
			word_list.append(w)

	total_score = c0_score + c2_score
	for k in dict1:
		v = blosum_matrix[k]
		if(k != c1_match and ((total_score + v) >= threshold)):
			if(word[1] == k[0]):
				k = k[1]
			else:
				k = k[0]
			w = word[0] + k + word[2]
			word_list.append(w)

	total_score = c0_score + c1_score
	for k in dict2:
		v = blosum_matrix[k]
		# print(k + " " + str(v))
		if(k != c2_match and ((total_score + v) >= threshold)):
			if(word[2] == k[0]):
				k = k[1]
			else:
				k = k[0]
			w = word[0] + word[1] + k
			word_list.append(w)

	return word_list


def construct_blosum_matrix():
	matrix_dict = {'AA':4, 'AC':0, 'AD':-2, 'AE':-1, 'AF':-2, 'AG':0, 'AH':-2, 'AI':-1, 'AK':-1, 'AL':-1, 'AM':-1, 'AN':-2, 'AP':-1, 'AQ':-1, 'AR':-1, 'AS':1,  
		'AT':0, 'AV':0, 'AW':-3, 'AY':-2, 'CC':9, 'CD':-3, 'CE':-4, 'CF':-2, 'CG':-3, 'CH':-3, 'CI':-1, 'CK':-3, 'CL':-1, 'CM':-1, 'CN':-3, 'CP':-3, 'CQ':-3,
		'CR':-3, 'CS':-1, 'CT':-1, 'CV':-1, 'CW':-2, 'CY':-2, 'DD':6,'DE':2, 'DF':-3, 'DG':-1, 'DH':-1, 'DI':-3, 'DK':-1, 'DL':-4, 'DM':-3, 'DN':1, 'DP':-1, 
		'DQ':0, 'DR':-2, 'DS':0, 'DT':-2, 'DV':-3, 'DW':-4, 'DY':-3, 'EE':5, 'EF':-3, 'EG':-2, 'EH':0, 'EI':-3, 'EK':1, 'EL':-3, 'EM':-2, 'EN':0, 'EP':-1, 
		'EQ':2, 'ER':0, 'ES':0, 'ET':-1, 'EV':-2, 'EW':-3, 'EY':-2, 'FF':6, 'FG':-3, 'FH':-1, 'FI':0, 'FK':-3, 'FL':0, 'FM':0, 'FN':-3, 'FP':-4, 'FQ':-3,
		'FR':-3, 'FS':-2, 'FT':-2, 'FV':-1, 'FW':1, 'FY':3, 'GG':6, 'GH':-2, 'GI':-4, 'GK':-2, 'GL':-4, 'GM':-3, 'GN':0, 'GP':-2, 'GQ':-2, 'GR':-2, 'GS':0, 
		'GT':-2, 'GV':-3, 'GW':-2, 'GY':-3, 'HH':8, 'HI':-3, 'HK':-1, 'HL':-3, 'HM':-2, 'HN':1, 'HP':-2, 'HQ':0, 'HR':0, 'HS':-1, 'HT':-2, 'HV':-3, 'HW':-2, 
		'HY':2, 'II':4, 'IK':-3, 'IL':2, 'IM':1, 'IN':-3, 'IP':-3, 'IQ':-3, 'IR':-3, 'IS':-1, 'IT':-1, 'IV':3, 'IW':-3, 'IY':-1, 'KK':5, 'KL':-2, 'KM':-1, 
		'KN':0, 'KP':-1, 'KQ':1, 'KR':2, 'KS':0, 'KT':-1, 'KV':-2, 'KW':-3, 'KY':-2, 'LL':4, 'LM':2, 'LN':-3, 'LP':-3, 'LQ':-2, 'LR':-2, 'LS':-2, 'LT':-1, 
		'LV':1, 'LW':-2, 'LY':-1, 'MM':5, 'MN':-2, 'MP':-2, 'MQ':0, 'MR':-1, 'MS':-1, 'MT':-1, 'MV':1, 'MW':-1, 'MY':-1, 'NN':6, 'NP':-2, 'NQ':0, 'NR':0, 
		'NS':1, 'NT':0, 'NV':-3, 'NW':-4, 'NY':-2, 'PP':7, 'PQ':-1, 'PR':-2, 'PS':-1, 'PT':-1, 'PV':-2, 'PW':-4, 'PY':-3, 'QQ':5, 'QR':1, 'QS':0, 'QT':-1, 
		'QV':-2, 'QW':-2, 'QY':-1, 'RR':5, 'RS':-2, 'RT':-2, 'RV':-3, 'RW':-3, 'RY':-2, 'SS':4, 'ST':1, 'SV':-2, 'SW':-3, 'SY':-2, 'TT':5, 'TV':0, 'TW':-2, 
		'TY':-2, 'VV':4, 'VW':-3, 'VY':-1, 'WW':11, 'WY':2, 'YY':7}

	return matrix_dict 

if __name__ == '__main__':
	args = sys.argv 
	file = args[1]
	threshold = args[2]
	threshold = int(threshold)

	query,subject = extract_params(file)
	words_q,words_s = split_query_subject(query, subject)
	
	matrix_dict = construct_blosum_matrix()
	dict_q, dict_s = record_indices(words_q, words_s)
	
	hsp = align(query, subject, words_q, words_s, dict_q, dict_s, matrix_dict, threshold)




