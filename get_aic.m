function [aic] = get_aic(negLL, n_params)
    % always change ModelDimen to the right dimensional space you are
    % testing
    ModelDimen = 1
    if ModelDimen == 1
      n_params = n_params - 2 % 1D    
    elseif ModelDimen == 2
        n_params = n_params - 4 % 2D
    end
    
    aic = 2*n_params + 2*negLL
end
