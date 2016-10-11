<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Nodal tetrahedron integration, electrostatic points, P = 6</description>
    <executable>NodalDemo</executable>
    <parameters>--order 6 --type 26 --integral</parameters>
    <metrics>
        <metric type="Linf" id="1">
            <value tolerance="1e-12">2.69479e-07</value>
        </metric>
    </metrics>
</test>
