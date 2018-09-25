<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh (v2.2) linear tri square, order 9 output</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list SquareTriLinear_2_2.msh SquareTri_2_2.xml:xml:test:order=9</parameters>
    <files>
        <file description="Input File">SquareTriLinear_2_2.msh</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
