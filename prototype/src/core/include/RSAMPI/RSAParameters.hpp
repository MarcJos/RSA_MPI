#pragma once

namespace rsa_mpi
{
  struct RSAParticleParameters
  {
    double radius;
    double volume_fraction;
    std::string type;
  };
}

namespace YAML
{
  using exaDEM::RSAParticleParameters;

  template <> struct convert<RSAParticleParameters>
  {
    static bool decode(const Node &node, RSAParticleParameters &v)
    {
      if (node.size() != 3)
      {
        return false;
      }
      v.radius = node[0].as<double>();
      v.volume_fraction = node[1].as<double>();
      v.type = node[2].as<std::string>();
      return true;
    }
  };
}
