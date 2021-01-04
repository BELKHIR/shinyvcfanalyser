awk '
    /def / {
            split($0,a,"[ |(]") 
            curModel = a[2]
            models[curModel]++
            }
    /"""/  {
        getline;
        #print curModel"\t"$0
        idparam = 0
        while( index($0, "\"\"\"") == 0)
        {   
            gsub(/=params/, "", $0);
            n = split($0, b, ":")  
            gsub(/ /, "", b[1]);
            if ((length(b)>= 2) && (b[1] != "n1,n2") && (b[1] != "pts"))
            { 
                idparam++; 
                gsub(/^ /, "", b[2]); gsub(/$ /, "", b[2]); 
                params[curModel][idparam] = b[1]"\t"b[2] 
                nbparams[curModel] = idparam 
            }
            getline
        }
    }     
    END {
        print "Model\tlib\thelp"
        for (model in models)
        {
            for(param=1; param<= nbparams[model]; param++ ) print model"\t"params[model][param]
        }
    }   
    '  ../momentsmodels/moments_models_2pop.py  > models_params.txt
