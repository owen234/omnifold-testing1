import numpy as np



def simple_pdf1a( theta0, theta_unknown_S, train_train_det_within_sphere, train_train_gen_within_sphere, iterations, scales, radius, verbose=0, compute_radii=True ):

    weights = np.empty(shape=(iterations, 2, len(theta0)))
    # shape = (iteration, step, event)
    push_weights_for_output = np.empty(shape=(iterations, len(theta0)))

    theta0_G = theta0[:,0]
    theta0_S = theta0[:,1]

    ngen_train = len(theta0)
    ngen_true  = len(theta_unknown_S)
    ndim = theta0_G.shape[1]

    labels0 = np.zeros(len(theta0))
    labels_unknown = np.ones(len(theta_unknown_S))
    labels_unknown_step2 = np.ones(len(theta0_G))


    if verbose :
        print("\n\n")
        print("  ======== simple_pdf1a\n\n")
        print("  shape of theta0_S : %s" % str(np.shape(theta0_S)) )
        print("  shape of theta0_G : %s" % str(np.shape(theta0_G)) )
        print("  shape of theta_unknown_S : %s" % str(np.shape(theta_unknown_S)) )
        print("\n  iterations = %d\n" % iterations )
        print("  scales :", scales )
        print("  radius = %.2f\n" % radius )
        print("  compute radii : ", compute_radii )
        print("  ndim : %d" % ndim )
        print("\n\n")

    # initial iterative weights are ones
    weights_pull = np.ones(len(theta0_S))
    weights_push = np.ones(len(theta0_S))



    #global train_train_det_within_sphere
    #global train_train_gen_within_sphere

    if compute_radii :

       if verbose :
          print('\n\n Calculating distance between all pairs of MC events')

          #-- create outside here
          #train_train_det_within_sphere = np.zeros( shape=(ngen_train ,ngen_train), dtype=bool )
          #train_train_gen_within_sphere = np.zeros( shape=(ngen_train ,ngen_train), dtype=bool )

          delta_det = np.zeros( shape=(ngen_train, ndim) )
          delta_gen = np.zeros( shape=(ngen_train, ndim) )


          for pi in range( ngen_train ) :
             if (pi%(ngen_train/10) == 0 ) : print('  %12d / %12d' % (pi, ngen_train) )
             dist2_det = np.zeros( shape=(ngen_train) )
             dist2_gen = np.zeros( shape=(ngen_train) )
             for di in range( ndim ) :
                 delta_det[:,di] = (theta0_S[:,di] - theta0_S[pi,di]) / scales[di]
                 delta_gen[:,di] = (theta0_G[:,di] - theta0_G[pi,di]) / scales[di]
                 dist2_det = dist2_det + delta_det[:,di] * delta_det[:,di]
                 dist2_gen = dist2_gen + delta_gen[:,di] * delta_gen[:,di]
             train_train_det_within_sphere[pi,:] = ( dist2_det < radius )
             train_train_gen_within_sphere[pi,:] = ( dist2_gen < radius )

          print('  --- Done calculating ')


    if verbose :
        print(' train_train_det_within_sphere ')
        print( train_train_det_within_sphere )



    if verbose :
        print('\n\n Calculating pdf from data, detector features')

    train_true_det_count_within_sphere = np.zeros( shape=(ngen_train), dtype=int )

    delta_det = np.zeros( shape=(ngen_true, ndim) )


    for pi in range( ngen_train ) :
        dist2_det = np.zeros( shape=(ngen_true) )
        for di in range( ndim ) :
            delta_det[:,di] = (theta_unknown_S[:,di] - theta0_S[pi,di]) / scales[di]
            dist2_det = dist2_det + delta_det[:,di] * delta_det[:,di]
        train_true_det_count_within_sphere[pi] = ( dist2_det < radius ).sum()

    true_det_simple_pdf_at_train_points = train_true_det_count_within_sphere / train_true_det_count_within_sphere.sum()

    if verbose :
       print('true_det_simple_pdf_at_train_points')
       print(true_det_simple_pdf_at_train_points)




    return_dict = {}




    for i in range(iterations):

        if (verbose>0):
            print("\nITERATION: {}\n".format(i + 1))
            pass

        # STEP 1: classify Sim. (which is reweighted by weights_push) to Data
        # weights reweighted Sim. --> Data

        if (verbose>0):
            print("   -- ITERATION %d  STEP 1\n" % (i+1) )
            pass

        print(" weights_push at the beginning")
        print( weights_push )



        train_det_simple_pdf_at_train_points = np.zeros( shape=(ngen_train) )

        for pi in range( ngen_train ) :

            train_det_simple_pdf_at_train_points[pi] = weights_push[ train_train_det_within_sphere[pi] ].sum()

        train_det_simple_pdf_at_train_points = train_det_simple_pdf_at_train_points / train_det_simple_pdf_at_train_points.sum()

        train_det_simple_pdf_at_train_points = np.clip( train_det_simple_pdf_at_train_points, 1e-11, 1e6)

        step1_output_weights = true_det_simple_pdf_at_train_points / train_det_simple_pdf_at_train_points

        weights_pull = weights_push * step1_output_weights

        weights[i, :1, :] = step1_output_weights



#       train_history_step1 = model.fit(X_train_1,
#                 Y_train_1,
#                 epochs=this_epochs,
#                 batch_size=batch_size_setval,
#                 validation_data=(X_test_1, Y_test_1),
#                 verbose=verbose)

#       step1_output_weights = reweight(theta0_S,model)

#       weights_pull = weights_push * step1_output_weights

#       weights[i, :1, :] = step1_output_weights

#       return_dict["train-hist-step1-iter%d" % i] = train_history_step1











        # STEP 2: classify Gen. to reweighted Gen. (which is reweighted by weights_pull)
        # weights Gen. --> reweighted Gen.

        if (verbose>0):
            print("\n   -- ITERATION %d  STEP 2\n" % (i+1) )
            pass



        train_gen_simple_pdf_at_train_points_push_weight = np.zeros( shape=(ngen_train) )

        for pi in range( ngen_train ) :

            train_gen_simple_pdf_at_train_points_push_weight[pi] = weights_push[ train_train_gen_within_sphere[pi] ].sum()

        train_gen_simple_pdf_at_train_points_push_weight = train_gen_simple_pdf_at_train_points_push_weight / train_gen_simple_pdf_at_train_points_push_weight.sum()





        train_gen_simple_pdf_at_train_points_pull_weight = np.zeros( shape=(ngen_train) )

        for pi in range( ngen_train ) :

            train_gen_simple_pdf_at_train_points_pull_weight[pi] = weights_pull[ train_train_gen_within_sphere[pi] ].sum()

        train_gen_simple_pdf_at_train_points_pull_weight = train_gen_simple_pdf_at_train_points_pull_weight / train_gen_simple_pdf_at_train_points_pull_weight.sum()




        train_gen_simple_pdf_at_train_points_push_weight = np.clip( train_gen_simple_pdf_at_train_points_push_weight, 1e-11, 1e6)

        step2_output_weights = train_gen_simple_pdf_at_train_points_pull_weight / train_gen_simple_pdf_at_train_points_push_weight

        weights_push = weights_push * step2_output_weights

        push_weights_for_output[i] = weights_push

        weights[i, 1:2, :] = step2_output_weights





#       step2_output_weights = reweight(theta0_G,model)

#       weights_push = step2_output_weights * weights_push

#       push_weights_for_output[i] = weights_push

#       weights[i, 1:2, :] = step2_output_weights

#       if save_step2_model :
#           model_output_dir = "%s/of-step2-iter%02d-model" % (output_dir, i)
#           print("\n +++ Saving step 2, iteration %d model in %s" % (i, model_output_dir) )
#           tf.keras.models.save_model( model, model_output_dir )


#       pass


    return_dict["weights"] = weights

    return_dict["push_weights"] = push_weights_for_output

    return_dict["final_push_weights"] = weights_push

    return return_dict





