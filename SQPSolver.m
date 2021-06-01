function [x,z,Hist] = SQPSolver(x0,obj,con,l,u,cl,cu, options)
    if nargin <8
        options = struct();
    end
    if sum(strcmp(fieldnames(options), 'log')) == 1
        if islogical(options.log)
            log = options.log;
        else
            error('The option log should be a boolean')
        end
    else
        log = 0;
    end
    if sum(strcmp(fieldnames(options), 'method')) == 1
        if ismember(options.method, {'SQP', 'SQP_ls', 'SQP_trust'})
            method = options.method;
        else
            error('The option method should be one of {SQP, SQP_ls, SQP_trust}')
        end
    else
        method = 'SQP_ls';
    end
    if sum(strcmp(fieldnames(options), 'infesibility_handling')) == 1
       if islogical(options.infesibility_handling)
           infesibility_handling = options.infesibility_handling;
       else
           error('The options infesibility_handling should be a boolean')
       end
    else
        infesibility_handling = true;
    end
    if sum(strcmp(fieldnames(options), 'precision')) == 1
        if isnumeric(options.precision) && floor(options.precision)==options.precision && options.precision>0 && options.precision<6
            precision = options.precision;
        else
            error('The option precision should be a positive integer below 6.')
        end
    else
        precision = 3;
    end
    if sum(strcmp(fieldnames(options), 'subsolver')) == 1
        if ismember(options.subsolver , {'own solver', 'quadprog'})
            if strcmp(options.subsolver,'own solver')
                %true is own solver
                subsolver = true;
            else
                subsolver = false;
            end
        else
            msg = ['The option subsolver should be either be ',  char(39),  'own solver', char(39) ,' or ', char(39), 'quadprog' ,char(39)];
            error(msg)
        end
    else
        subsolver = true;
    end
    
    
    
    if method == "SQP" && infesibility_handling==0
        [x,z,Hist] = SQP(x0,obj,con,l,u,cl,cu,log,subsolver,precision);
        
    elseif method == "SQP_ls" && infesibility_handling==0
        [x,z,Hist] = SQP_ls(x0,obj,con,l,u,cl,cu,log,subsolver,precision);
        
    elseif method == "SQP" && infesibility_handling==1
        [x,z,Hist] = SQP_infes(x0,obj,con,l,u,cl,cu,log,precision);
        
    elseif method == "SQP_ls" && infesibility_handling==1
        [x,z,Hist] = SQP_ls_infes(x0,obj,con,l,u,cl,cu,log,precision);
        
    elseif method == "SQP_trust"
        [x,z,Hist] = SQP_trust(x0,obj,con,l,u,cl,cu,log,precision);
    
    end

end

