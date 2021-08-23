function ctrl_state = state2ctrl_state(X,ctrl_state)
    ctrl_state(3) = X(1);
    ctrl_state(6) = X(2);
    ctrl_state(7:9) = [X(3); X(5); X(7)];
    ctrl_state(10:12) = [X(4); X(6); X(8)];
    ctrl_state(13:16) = [X(9); X(10); X(11); X(12)];
end