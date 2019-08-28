function [s] = moresorensen( g, H, Delta, eps_D )
  
   theta    = 0.0001;
   n        = length(g);
   value    = 0;
   gnorm    = norm(g);
   g2D      = gnorm/Delta;
   Hnorminf = norm(H, inf );
   HnormF   = norm(H, 'fro' );
   lower    = max( 0, g2D - min( Hnorminf, HnormF ) );
   upper    = max( 0, g2D + min( Hnorminf, HnormF ) );
   Dlower   = ( 1 - eps_D ) * Delta;
   Dupper   = ( 1 + eps_D ) * Delta;

   if ( lower < 1.0e-12 )
      lambda = 0;
   else
      lambda = max( sqrt( lower * upper ), lower + theta * ( upper - lower) );
   end
   
   for i = 1:100
      new_lambda = -1;
      [ R, p ] = chol( H + lambda * speye( n ) );
      if ( p == 0 )
         s     = - R \ ( R' \ g );
         norms = norm( s );
         if ( ( lambda < 1.0e-12 & norms <= Dupper ) || ( norms >= Dlower & norms<= Dupper ) )
            break;
         end
         w          = R' \ s;
         normw2     = norm( w )^2;
         new_lambda = lambda + ( ( norms - Delta ) / Delta ) * ( norms^2 / normw2 );
         if ( norms > Dupper )
            lower = lambda;
         else
            upper = lambda;
         end
         theta_range = theta * ( upper - lower );
         if ( new_lambda > lower + theta_range & new_lambda < upper - theta_range )
            lambda = new_lambda;
         else
            lambda = max( sqrt( lower * upper ), lower + theta_range );
         end
      else
         lower  = lambda;
         lambda = max( sqrt( lower * upper ), lower + theta * ( upper - lower) );
      end
   end

   return