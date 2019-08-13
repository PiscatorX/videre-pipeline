

with open("ko_map.tsv") as fp:
    for line in fp:
        #k141_1042899', 'K01084', 'G6PC\tglucose-6-phosphatase \tEC:3.1.3.9\n']
        contig, K_number, annotation =  line.split("\t",2)
        definition = annotation.strip().split("\t")[:2]
        identifier, name = definition[:2]
        ec_number  = ''
        if name.endswith("]"):
            name, ec_number = name.rstrip("]").split("[EC:")
        elif len(definition) == 3:
            ec_number = annotation[2]
        print(K_number,  name,  ec_number)
        
