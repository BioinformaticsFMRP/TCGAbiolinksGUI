$(function() {
  $(document).on({

    'shiny:connected': function(event) {
      $('form.well').fadeOut(3000).fadeIn(2000);
    },

    'shiny:disconnected': function(event) {
      //alert('Disconnected! The web socket state is ' + event.socket.readyState);
      close();
    },

    'shiny:busy': function(event) {
       console.log('Busy ' + new Date());
      //$('#busyModal').modal('show');
    },

    'shiny:idle': function(event) {
      console.log('Idle ' + new Date());
      //$('#busyModal').modal('hide');
    },

    'shiny:recalculating': function(event) {
      console.log('An output is being recalculated... ' + new Date());
      $('#busyModal').modal('show');
    },

    'shiny:recalculated': function(event) {
      console.log('An output has been recalculated! ' + new Date());
      $('#busyModal').modal('hide');
    }

  });


  Shiny.addCustomMessageHandler('special', function(message) {
    //
  });
});
