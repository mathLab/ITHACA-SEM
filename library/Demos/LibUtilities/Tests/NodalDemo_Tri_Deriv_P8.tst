<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Nodal triangle derivative, electrostatic points, P = 8</description>
    <executable>NodalDemo</executable>
    <parameters>--order 8 --type 22 --deriv</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.000192327</value>
            <value variable="v" tolerance="1e-12">0.000192208</value>
        </metric>
    </metrics>
</test>
