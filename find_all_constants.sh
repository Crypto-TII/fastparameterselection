#!/bin/bash
# python3 src/estimate.py --fit --param "lambda" --attack 'usvp' --secret "binary" --error "3.19" --simpl 0
# python3 src/estimate.py --fit --param "lambda" --attack 'usvp' --secret "ternary" --error "3.19" --simpl 0
# python3 src/estimate.py --fit --param "lambda" --attack 'usvp' --secret "binary" --error "3.19" --simpl 1
# python3 src/estimate.py --fit --param "lambda" --attack 'usvp' --secret "ternary" --error "3.19" --simpl 1

# python3 src/estimate.py --fit --param "lambda" --attack 'bdd' --secret "binary" --error "3.19" --simpl 0
# python3 src/estimate.py --fit --param "lambda" --attack 'bdd' --secret "ternary" --error "3.19" --simpl 0
# python3 src/estimate.py --fit --param "lambda" --attack 'bdd' --secret "binary" --error "3.19" --simpl 1
# python3 src/estimate.py --fit --param "lambda" --attack 'bdd' --secret "ternary" --error "3.19" --simpl 1

# python3 src/estimate.py --fit --param "n" --attack 'usvp' --secret "binary" --error "3.19" --simpl 0
# python3 src/estimate.py --fit --param "n" --attack 'usvp' --secret "ternary" --error "3.19" --simpl 0
# python3 src/estimate.py --fit --param "n" --attack 'usvp' --secret "binary" --error "3.19" --simpl 1
# python3 src/estimate.py --fit --param "n" --attack 'usvp' --secret "ternary" --error "3.19" --simpl 1

python3 src/estimate.py --fit --param "n" --attack 'bdd' --secret "binary" --error "3.19" --simpl 0
python3 src/estimate.py --fit --param "n" --attack 'bdd' --secret "ternary" --error "3.19" --simpl 0
python3 src/estimate.py --fit --param "n" --attack 'bdd' --secret "binary" --error "3.19" --simpl 1
python3 src/estimate.py --fit --param "n" --attack 'bdd' --secret "ternary" --error "3.19" --simpl 1