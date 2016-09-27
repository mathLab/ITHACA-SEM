<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Nodal tetrahedron derivative, electrostatic points, P = 8</description>
    <executable>NodalDemo</executable>
    <parameters>--order 8 --type 26 --deriv</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.000519629</value>
            <value variable="v" tolerance="1e-12">0.000577377</value>
            <value variable="w" tolerance="1e-12">0.000519629</value>
        </metric>
    </metrics>
</test>
