% file MASS/man/cars93.Rd
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{cars93}
\alias{cars93}
\title{
  Data from 82 Cars on Sale in the USA in 1993
}
\description{
  The \code{cars93} data frame has 82 rows and 16 columns.
}
\usage{
cars93
}
\format{
  This data frame contains the following columns:
  \describe{

    \item{\code{MPG.city}}{
      City MPG (miles per US gallon by EPA rating).
    }
    \item{\code{MPG.highway}}{
      Highway MPG.
    }
    \item{\code{AirBags}}{
      Air Bags standard. Factor: none, driver only, or driver & passenger.
    }
    \item{\code{EngineSize}}{
      Engine size (litres).
    }
    \item{\code{Horsepower}}{
      Horsepower (maximum).
    }
    \item{\code{RPM}}{
      RPM (revs per minute at maximum horsepower).
    }
    \item{\code{Rev.per.mile}}{
      Engine revolutions per mile (in highest gear).
    }
    \item{\code{Fuel.tank.capacity}}{
      Fuel tank capacity (US gallons).
    }
    \item{\code{Passengers}}{
      Passenger capacity (persons)
    }
    \item{\code{Length}}{
      Length (inches).
    }
    \item{\code{Width}}{
      Width (inches).
    }
    \item{\code{Turn.circle}}{
      U-turn space (feet).
    }
    \item{\code{Rear.seat.room}}{
      Rear seat room (inches) (missing for 2-seater vehicles).
    }
    \item{\code{Luggage.room}}{
      Luggage capacity (cubic feet) (missing for vans).
    }
    \item{\code{Make}}{
      Combination of Manufacturer and Model (character).
    }

    \item{\code{price}}{
      Midrange Price (in \$1,000): average of \code{Min.Price} and
      \code{Max.Price}.
    }

  }
}
\details{
  Cars were selected at random from among 1993 passenger car models that
  were listed in both the \emph{Consumer Reports} issue and the
  \emph{PACE Buying Guide}.  Pickup trucks and Sport/Utility vehicles were
  eliminated due to incomplete information in the \emph{Consumer Reports}
  source.  Duplicate models (e.g., Dodge Shadow and Plymouth Sundance)
  were listed at most once.

  Further description can be found in Lock (1993).
}
\source{
  Lock, R. H. (1993)
  1993 New Car Data.
  \emph{Journal of Statistics Education}
  \bold{1}(1).
  \url{http://www.amstat.org/publications/jse/v1n1/datasets.lock.html}.  
}
\references{
  Venables, W. N. and Ripley, B. D. (1999)
  \emph{Modern Applied Statistics with S-PLUS.} Third
  Edition. Springer.
}
\keyword{datasets}
